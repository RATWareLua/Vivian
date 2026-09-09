--!native
-- chromosome.lua  (phase 2 -- the chromosome layer)
--
-- A chromosome is a packed strand built on genes (genome.lua):
--   [TELOMERE_L][CENTROMERE][gene 0]..[gene n-1][TELOMERE_R]
--
-- Telomeres: TTAGGG-style repeats (digits 3,3,0,2,2,2) capping both
-- ends; damage there is detectable (pattern mismatch) and repairable.
-- Centromere: marker + codon-coded header (chr id, flags, gene count,
-- generation counter) protected by its own Chaskey-12 tag. The generation
-- counter is the telomere-analog: it grows with every replication and
-- drives senescence / self-renewal (phase 4).
--
-- Genes are found by the gene scanner inside the interior region, so a
-- chromosome does not need to know gene offsets -- data segments stay
-- self-delimiting, like on a real chromosome.

local genome = require("./genome")

local byte, char, concat, sub, find =
	string.byte, string.char, table.concat, string.sub, string.find
local band, rshift =
	bit32.band, bit32.rshift

-- telomere repeat chunk: 2x TTAGGG = 12 bases = 3 bytes
local TELUNIT = "\242\227\14"
local CMARK = "\109\109\109"  -- 12 bases CTGA x3 (digits 1,2,3,1)
local CENBYTES = 18           -- marker (3) + 20 codons (15)

local function build_centromere(chr_id, flags, ngenes, generation)
	local raw = char(band(chr_id, 255), band(flags, 255),
		band(rshift(ngenes, 8), 255), band(ngenes, 255),
		band(rshift(generation, 8), 255), band(generation, 255))
	local c = genome.chaskey32(raw)
	local block = raw .. char(band(rshift(c, 24), 255), band(rshift(c, 16), 255),
		band(rshift(c, 8), 255), band(c, 255))
	return CMARK .. genome.pack_codons(genome.values_from_bytes(block))
end

local function encode(chr_id, data, opts)
	opts = opts or {}
	if type(chr_id) ~= "number" or chr_id ~= math.floor(chr_id)
		or chr_id < 0 or chr_id > 255 then
		return nil, "invalid chr_id (byte 0..255)"
	end
	if type(data) ~= "string" then return nil, "expected string" end
	local flags = opts.flags or 0
	if type(flags) ~= "number" or flags ~= math.floor(flags)
		or flags < 0 or flags > 255 then
		return nil, "invalid flags (byte 0..255)"
	end
	local units = opts.units or 4
	if type(units) ~= "number" or units ~= math.floor(units) or units < 1 then
		return nil, "invalid units (integer >= 1)"
	end
	local gene_raw = opts.gene_raw or 1024
	if type(gene_raw) ~= "number" or gene_raw ~= math.floor(gene_raw) or gene_raw < 16 then
		return nil, "invalid gene_raw (integer >= 16)"
	end
	local generation = opts.generation or 0
	if type(generation) ~= "number" or generation ~= math.floor(generation)
		or generation < 0 or generation > 65535 then
		return nil, "invalid generation (0..65535)"
	end
	local telo = TELUNIT:rep(units)
	local genes, ngenes = {}, 0
	for i = 1, #data, gene_raw do
		local g = genome.gene_encode(sub(data, i, i + gene_raw - 1),
			{ id = ngenes, mode = opts.mode or "dense", h = opts.h })
		if not g then return nil, "gene encode failed" end
		ngenes = ngenes + 1
		genes[ngenes] = g
	end
	local cen = build_centromere(chr_id, flags, ngenes, generation)
	return telo .. cen .. concat(genes) .. telo
end

-- parse returns a chromosome record; gene offsets are absolute (into
-- the strand), so repairs can splice in place
-- verify the centromere block starting at byte p (marker included) by
-- its CRC; used to validate auto-detected telomere size even when the
-- telomere itself is damaged
local function cen_valid_at(strand, p)
	local cenb = sub(strand, p, p + CENBYTES - 1)
	if #cenb < CENBYTES or sub(cenb, 1, 3) ~= CMARK then return false end
	local v = genome.unpack_codons(sub(cenb, 4), 20)
	local raw = genome.bytes_from_values(v, 6)
	local expected = (v[13] * 16 + v[14]) * 0x1000000
		+ (v[15] * 16 + v[16]) * 0x10000
		+ (v[17] * 16 + v[18]) * 0x100
		+ (v[19] * 16 + v[20])
	return genome.chaskey32(raw) == expected
end

local function parse(strand, opts)
	if type(strand) ~= "string" then return nil, "expected string" end
	-- auto-detect the telomere size from the centromere marker position:
	-- the marker is at a fixed offset (damage never shifts the layout),
	-- and the candidate is validated by the centromere CRC at that
	-- offset, so it works even when the telomere itself is damaged and
	-- phantom markers in payloads cannot pass the CRC check
	local p_cmark = find(strand, CMARK, 1, true)
	local units = (opts and opts.units) or 4
	if p_cmark and p_cmark > 1 and (p_cmark - 1) % #TELUNIT == 0
		and cen_valid_at(strand, p_cmark) then
		units = (p_cmark - 1) / #TELUNIT
	end
	local telo = TELUNIT:rep(units)
	local tb = #telo
	if #strand < 2 * tb + CENBYTES then return nil, "chromosome too short" end
	local telo_ok = sub(strand, 1, tb) == telo
		and sub(strand, #strand - tb + 1) == telo
	local cen_ok = false
	local id, flags, ngenes, generation
	local cenb = sub(strand, tb + 1, tb + CENBYTES)
	if sub(cenb, 1, 3) == CMARK then
		local v = genome.unpack_codons(sub(cenb, 4), 20)
		local raw = genome.bytes_from_values(v, 6)
		local expected = (v[13] * 16 + v[14]) * 0x1000000
			+ (v[15] * 16 + v[16]) * 0x10000
			+ (v[17] * 16 + v[18]) * 0x100
			+ (v[19] * 16 + v[20])
		if genome.chaskey32(raw) == expected then
			cen_ok = true
			id = byte(raw, 1)
			flags = byte(raw, 2)
			ngenes = byte(raw, 3) * 256 + byte(raw, 4)
			generation = byte(raw, 5) * 256 + byte(raw, 6)
		end
	end
	local genes = {}
	if cen_ok then
		local interior = sub(strand, tb + CENBYTES + 1, #strand - tb)
		local base = tb + CENBYTES
		local gs = genome.scan(interior)
		for i = 1, #gs do
			local g = gs[i]
			genes[i] = {
				id = g.id, mode = g.mode, rawlen = g.rawlen, packedlen = g.packedlen,
				offset = base + g.offset, size = g.size,
				data = g.data, crc_ok = g.crc_ok,
			}
		end
	end
	return {
		strand = strand,
		id = id, flags = flags, ngenes = ngenes, generation = generation,
		telomere_ok = telo_ok, cen_ok = cen_ok,
		genes = genes, telomere_bytes = tb, units = units,
	}
end

local function read(rec)
	if not rec or not rec.cen_ok then return nil, "centromere damaged" end
	local n = rec.ngenes
	if #rec.genes ~= n then return nil, "gene count mismatch" end
	local seen = {}
	for i = 1, n do
		local g = rec.genes[i]
		if seen[g.id] then return nil, "duplicate gene id" end
		seen[g.id] = g
	end
	local parts = {}
	for i = 0, n - 1 do
		local g = seen[i]
		if not g then return nil, "gene " .. i .. " missing" end
		if not g.crc_ok then return nil, "gene " .. i .. " damaged" end
		parts[i + 1] = g.data
	end
	return concat(parts)
end

local function set_generation(strand, generation, opts)
	local rec = parse(strand, opts)
	if not rec or not rec.cen_ok then return nil, "centromere damaged" end
	if type(generation) ~= "number" or generation ~= math.floor(generation)
		or generation < 0 or generation > 65535 then
		return nil, "invalid generation (0..65535)"
	end
	local tb = rec.telomere_bytes
	local cen = build_centromere(rec.id, rec.flags, rec.ngenes, generation)
	return sub(strand, 1, tb) .. cen .. sub(strand, tb + CENBYTES + 1)
end

return {
	encode = encode,
	parse = parse,
	read = read,
	set_generation = set_generation,
	build_centromere = build_centromere,
	TELUNIT = TELUNIT,
	CMARK = CMARK,
}

