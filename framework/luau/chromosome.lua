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
--
-- VIV14N: chr_opts.parity = m appends m Reed-Solomon parity genes after
-- the data genes and switches the centromere to revision 2 (raw block
-- id, flags, ngenes, generation, parity, reserved, rawlen). Any n of the
-- n+m genes reconstruct the payload; gene type bit 1 marks parity genes.

local genome = require("./genome")
local parity = require("./parity")

local byte, char, concat, sub, find =
	string.byte, string.char, table.concat, string.sub, string.find
local floor, min = math.floor, math.min
local band, rshift = bit32.band, bit32.rshift

-- telomere repeat chunk: 2x TTAGGG = 12 bases = 3 bytes
local TELUNIT = "\242\227\14"
local CMARK = "\109\109\109"  -- 12 bases CTGA x3 (digits 1,2,3,1)
local CENBYTES = 18           -- marker (3) + 20 codons (15)
local CEN2BYTES = 27          -- marker (3) + 32 codons (24)

local function build_centromere(chr_id, flags, ngenes, generation)
	local raw = char(band(chr_id, 255), band(flags, 255),
		band(rshift(ngenes, 8), 255), band(ngenes, 255),
		band(rshift(generation, 8), 255), band(generation, 255))
	local c = genome.chaskey32(raw)
	local block = raw .. char(band(rshift(c, 24), 255), band(rshift(c, 16), 255),
		band(rshift(c, 8), 255), band(c, 255))
	return CMARK .. genome.pack_codons(genome.values_from_bytes(block))
end

local function build_centromere2(chr_id, flags, ngenes, generation, np, rawlen)
	local raw = char(band(chr_id, 255), band(flags, 255),
		band(rshift(ngenes, 8), 255), band(ngenes, 255),
		band(rshift(generation, 8), 255), band(generation, 255),
		band(np, 255), 0,
		band(rshift(rawlen, 24), 255), band(rshift(rawlen, 16), 255),
		band(rshift(rawlen, 8), 255), band(rawlen, 255))
	local c = genome.chaskey32(raw)
	local block = raw .. char(band(rshift(c, 24), 255), band(rshift(c, 16), 255),
		band(rshift(c, 8), 255), band(c, 255))
	return CMARK .. genome.pack_codons(genome.values_from_bytes(block))
end

local function encode(chr_id, data, opts)
	opts = opts or {}
	if type(chr_id) ~= "number" or chr_id ~= floor(chr_id)
		or chr_id < 0 or chr_id > 255 then
		return nil, "invalid chr_id (byte 0..255)"
	end
	if type(data) ~= "string" then return nil, "expected string" end
	local flags = opts.flags or 0
	if type(flags) ~= "number" or flags ~= floor(flags)
		or flags < 0 or flags > 255 then
		return nil, "invalid flags (byte 0..255)"
	end
	local units = opts.units or 4
	if type(units) ~= "number" or units ~= floor(units) or units < 1 then
		return nil, "invalid units (integer >= 1)"
	end
	local gene_raw = opts.gene_raw or 1024
	if type(gene_raw) ~= "number" or gene_raw ~= floor(gene_raw) or gene_raw < 16 then
		return nil, "invalid gene_raw (integer >= 16)"
	end
	local generation = opts.generation or 0
	if type(generation) ~= "number" or generation ~= floor(generation)
		or generation < 0 or generation > 65535 then
		return nil, "invalid generation (0..65535)"
	end
	local np = opts.parity or 0
	if type(np) ~= "number" or np ~= floor(np) or np < 0 or np > 16 then
		return nil, "invalid parity (integer 0..16)"
	end
	local telo = TELUNIT:rep(units)
	local genes, ngenes = {}, 0
	for i = 1, #data, gene_raw do
		if ngenes >= 256 then   -- gene ids are one byte: 0..255
			return nil, "too many genes (max 256)"
		end
		local g, err = genome.gene_encode(sub(data, i, i + gene_raw - 1),
			{ id = ngenes, mode = opts.mode or "dense", h = opts.h })
		if not g then return nil, "gene encode failed" end
		ngenes = ngenes + 1
		genes[ngenes] = g
	end
	-- generation is stamped by set_generation / replication, encode writes 0
	local cen = build_centromere(chr_id, flags, ngenes, 0)
	if np == 0 then
		return telo .. cen .. concat(genes) .. telo
	end
	if ngenes == 0 then return nil, "parity needs a non-empty payload" end
	if ngenes + np > 255 then return nil, "too many genes with parity (max 255)" end

	-- data shards zero-padded to gene_raw -> systematic RS parity shards
	local shards_in = {}
	for i = 1, ngenes do
		local off = (i - 1) * gene_raw
		local chunk = sub(data, off + 1, off + gene_raw)
		shards_in[i] = chunk .. string.rep("\0", gene_raw - #chunk)
	end
	local shards, perr = parity.encode(shards_in, np, gene_raw)
	if not shards then return nil, perr end
	local par_genes = {}
	for j = 1, np do
		local g, perr2 = genome.gene_encode(shards[ngenes + j],
			{ id = ngenes + j - 1, type = 2, mode = opts.mode or "dense", h = opts.h })
		if not g then return nil, "parity gene encode failed: " .. tostring(perr2) end
		par_genes[j] = g
	end
	local cen2 = build_centromere2(chr_id, flags, ngenes + np, 0, np, #data)
	return telo .. cen2 .. concat(genes) .. concat(par_genes) .. telo
end

-- parse returns a chromosome record; gene offsets are absolute (into
-- the strand), so repairs can splice in place
-- parse either centromere revision at position p (1-based, marker
-- included); returns the fields plus version (2 = VIV14N, 1 = VIV1),
-- or nil when the block is not a valid tagged centromere
local function cen_parse(strand, p)
	if p + CENBYTES - 1 > #strand or sub(strand, p, p + 2) ~= CMARK then
		return nil
	end
	if p + CEN2BYTES - 1 <= #strand then
		local v = genome.unpack_codons(sub(strand, p + 3, p + 2 + 24), 32)
		local raw = genome.bytes_from_values(v, 12)
		local expected = (v[25] * 16 + v[26]) * 0x1000000
			+ (v[27] * 16 + v[28]) * 0x10000
			+ (v[29] * 16 + v[30]) * 0x100
			+ (v[31] * 16 + v[32])
		if genome.chaskey32(raw) == expected then
			return {
				version = 2,
				id = byte(raw, 1), flags = byte(raw, 2),
				ngenes = byte(raw, 3) * 256 + byte(raw, 4),
				generation = byte(raw, 5) * 256 + byte(raw, 6),
				parity = byte(raw, 7),
				rawlen = byte(raw, 9) * 0x1000000 + byte(raw, 10) * 0x10000
					+ byte(raw, 11) * 0x100 + byte(raw, 12),
			}
		end
	end
	local v = genome.unpack_codons(sub(strand, p + 3, p + CENBYTES - 1), 20)
	local raw = genome.bytes_from_values(v, 6)
	local expected = (v[13] * 16 + v[14]) * 0x1000000
		+ (v[15] * 16 + v[16]) * 0x10000
		+ (v[17] * 16 + v[18]) * 0x100
		+ (v[19] * 16 + v[20])
	if genome.chaskey32(raw) == expected then
		return {
			version = 1,
			id = byte(raw, 1), flags = byte(raw, 2),
			ngenes = byte(raw, 3) * 256 + byte(raw, 4),
			generation = byte(raw, 5) * 256 + byte(raw, 6),
			parity = 0, rawlen = 0,
		}
	end
	return nil
end

local function parse(strand, opts)
	if type(strand) ~= "string" then return nil, "expected string" end
	-- auto-detect the telomere size from the centromere marker position:
	-- the marker is at a fixed offset (damage never shifts the layout),
	-- and the candidate is validated by the centromere tag at that
	-- offset, so it works even when the telomere itself is damaged and
	-- phantom markers in payloads cannot pass the tag check
	local p_cmark = find(strand, CMARK, 1, true)
	local units = (opts and opts.units) or 4
	if p_cmark and p_cmark > 1 and (p_cmark - 1) % #TELUNIT == 0
		and cen_parse(strand, p_cmark) then
		units = (p_cmark - 1) / #TELUNIT
	end
	local telo = TELUNIT:rep(units)
	local tb = #telo
	if #strand < 2 * tb + CENBYTES then return nil, "chromosome too short" end
	local telo_ok = sub(strand, 1, tb) == telo
		and sub(strand, #strand - tb + 1) == telo
	local f = cen_parse(strand, tb + 1)
	if f and f.version == 2 and #strand < 2 * tb + CEN2BYTES then
		return nil, "chromosome too short"
	end
	local genes = {}
	if f then
		local cen_bytes = (f.version == 2) and CEN2BYTES or CENBYTES
		local interior = sub(strand, tb + cen_bytes + 1, #strand - tb)
		local base = tb + cen_bytes
		local gs = genome.scan(interior)
		for i = 1, #gs do
			local g = gs[i]
			genes[i] = {
				id = g.id, type = g.type, mode = g.mode, rawlen = g.rawlen,
				packedlen = g.packedlen,
				offset = base + g.offset, size = g.size,
				data = g.data, crc_ok = g.crc_ok,
			}
		end
	end
	return {
		strand = strand,
		id = f and f.id or nil, flags = f and f.flags or nil,
		ngenes = f and f.ngenes or nil, generation = f and f.generation or nil,
		parity = f and f.parity or 0, rawlen = f and f.rawlen or 0,
		cen_version = f and f.version or 0,
		telomere_ok = telo_ok, cen_ok = f ~= nil,
		genes = genes, telomere_bytes = tb, units = units,
	}
end

local function read(rec)
	if not rec or not rec.cen_ok then return nil, "centromere damaged" end
	local m = rec.parity or 0

	if m == 0 then
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

	-- VIV14N: reconstruct data shards from whatever genes survived
	local total = rec.ngenes
	local rawlen = rec.rawlen
	if total > 255 or m > total or rawlen == 0 then
		return nil, "invalid parity layout"
	end
	local n = total - m
	local seen = {}
	for i = 1, #rec.genes do
		local g = rec.genes[i]
		if g.id < 0 or g.id >= total then return nil, "gene id out of range" end
		if seen[g.id] then return nil, "duplicate gene id" end
		seen[g.id] = g
	end
	local all_data = true
	for k = 0, n - 1 do
		local g = seen[k]
		if not g or not g.crc_ok then all_data = false break end
	end
	if all_data then
		local parts, pos = {}, 0
		for k = 0, n - 1 do
			if pos >= rawlen then break end
			local g = seen[k]
			local take = min(#g.data, rawlen - pos)
			parts[#parts + 1] = sub(g.data, 1, take)
			pos = pos + take
		end
		local out = concat(parts)
		if #out ~= rawlen then return nil, "payload length mismatch" end
		return out
	end

	-- an erasure must be covered by parity: the shard length comes from
	-- a readable parity gene
	local shard_len = 0
	for k = n, total - 1 do
		local g = seen[k]
		if g and g.crc_ok then shard_len = g.rawlen break end
	end
	if shard_len == 0 or rawlen > n * shard_len then return nil, "too few genes" end
	local refs, present = {}, {}
	for k = 0, total - 1 do
		local g = seen[k]
		if g and g.crc_ok then
			local take = min(#g.data, shard_len)
			refs[k + 1] = sub(g.data, 1, take) .. string.rep("\0", shard_len - take)
			present[k + 1] = true
		end
	end
	local decoded, err = parity.decode(refs, present, n, m, shard_len)
	if not decoded then return nil, err end
	local parts, pos = {}, 0
	for k = 0, n - 1 do
		if pos >= rawlen then break end
		local take = min(shard_len, rawlen - pos)
		parts[#parts + 1] = sub(decoded[k + 1], 1, take)
		pos = pos + take
	end
	local out = concat(parts)
	if #out ~= rawlen then return nil, "payload length mismatch" end
	return out
end

local function set_generation(strand, generation, opts)
	local rec = parse(strand, opts)
	if not rec or not rec.cen_ok then return nil, "centromere damaged" end
	if type(generation) ~= "number" or generation ~= floor(generation)
		or generation < 0 or generation > 65535 then
		return nil, "invalid generation (0..65535)"
	end
	local tb = rec.telomere_bytes
	local cen, cen_bytes
	if rec.cen_version == 2 then
		cen = build_centromere2(rec.id, rec.flags, rec.ngenes, generation,
			rec.parity, rec.rawlen)
		cen_bytes = CEN2BYTES
	else
		cen = build_centromere(rec.id, rec.flags, rec.ngenes, generation)
		cen_bytes = CENBYTES
	end
	return sub(strand, 1, tb) .. cen .. sub(strand, tb + cen_bytes + 1)
end

return {
	encode = encode,
	parse = parse,
	read = read,
	set_generation = set_generation,
	build_centromere = build_centromere,
	build_centromere2 = build_centromere2,
	TELUNIT = TELUNIT,
	CMARK = CMARK,
	CENBYTES = CENBYTES,
	CEN2BYTES = CEN2BYTES,
}
