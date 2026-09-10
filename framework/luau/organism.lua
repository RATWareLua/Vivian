--!native
-- organism.lua  (phase 5 -- the organism: a diploid genome of N chromosomes)
--
-- An organism is a cell whose homologs are WHOLE GENOMES: tables of
-- chromosome strands, one homologous pair per chromosome. All repair,
-- replication and self-renewal machinery works per chromosome through
-- the cell.lua primitives.
--
-- Data-level evolution API (distinct from bit-level damage):
--   * mutate(): flips random BITS in the RAW data of a random
--     chromosome and re-encodes it -- the mutation is "fixed" into the
--     genome with a fresh Chaskey tag, so the organism can read its
--     own mutated genes afterwards (bit-level damage() would instead
--     be fought off by the repair machinery)
--   * cross(parentA, parentB): sexual reproduction -- each child
--     homolog is a per-gene mosaic of the parents' homologs (gene
--     level crossing over), re-encoded with fresh tags
--
-- Serialization: a versioned container ("VIV1") holding the canonical
-- genome (homolog 1) with per-chromosome coding options, so a saved
-- organism can be reloaded and keeps its layout for mutate/cross.
--
-- Requires bit32 (Luau, Lua 5.2+), cell, chromosome, genome, dna.

local cell = require("./cell")
local chromosome = require("./chromosome")

local byte, char, concat, sub =
	string.byte, string.char, table.concat, string.sub
local floor = math.floor

local band, bxor, lshift, rshift, bor =
	bit32.band, bit32.bxor, bit32.lshift, bit32.rshift, bit32.bor
local unpack = table.unpack or unpack

local MAGIC = "VIV1"
local MAGIC14 = "VIV14"
local MAGIC14N = "VIV14N"
local MAGIC_BANSHEE = "VIV14NB4NSH33"

local function bytes_to_string(t)
	local res = {}
	for i = 1, #t, 2048 do
		local j = (i + 2047 <= #t) and (i + 2047) or #t
		res[#res + 1] = char(unpack(t, i, j))
	end
	return concat(res)
end

local function opts_equal(o1, o2)
	o1, o2 = o1 or {}, o2 or {}
	return (o1.gene_raw or 1024) == (o2.gene_raw or 1024)
		and (o1.mode or "dense") == (o2.mode or "dense")
		and (o1.units or 4) == (o2.units or 4)
		and (o1.h or 3) == (o2.h or 3)
		and (o1.flags or 0) == (o2.flags or 0)
end

local function new(specs, opts)
	opts = opts or {}
	if type(specs) ~= "table" or #specs == 0 then
		return nil, "expected non-empty specs array"
	end
	if #specs > 255 then return nil, "too many chromosomes (max 255)" end
	local ids, chr_opts, gA, gB = {}, {}, {}, {}
	for i = 1, #specs do
		local sp = specs[i]
		if type(sp) ~= "table" then return nil, "spec " .. i .. " is not a table" end
		local cid = sp.id or i - 1
		if type(cid) ~= "number" or cid ~= floor(cid) or cid < 0 or cid > 255 then
			return nil, "invalid chr_id in spec " .. i
		end
		for j = 1, i - 1 do
			if ids[j] == cid then return nil, "duplicate chr_id " .. cid end
		end
		local strand, err = chromosome.encode(cid, sp.data or "", sp.opts or {})
		if not strand then return nil, "chromosome " .. i .. ": " .. tostring(err) end
		ids[i] = cid
		chr_opts[i] = sp.opts or {}
		gA[i], gB[i] = strand, strand
	end
	return {
		chr_ids = ids,
		chr_opts = chr_opts,
		homologs = { gA, gB },
		generation = 0,
		max_gen = opts.max_gen or 60,
		stem = false,
		dead = false,
		stem_source = nil,
	}
end

local function stem(specs, opts)
	local o, err = new(specs, opts)
	if not o then return nil, err end
	o.stem = true
	return o
end

local function attach_stem(o, stemOrg)
	o.stem_source = stemOrg
	return o
end

local function renew(o, stemOrg)
	if not stemOrg or not stemOrg.homologs then return nil, "no stem organism" end
	o.homologs = {
		{ unpack(stemOrg.homologs[1], 1, #stemOrg.chr_ids) },
		{ unpack(stemOrg.homologs[2], 1, #stemOrg.chr_ids) },
	}
	o.chr_ids = stemOrg.chr_ids
	o.chr_opts = stemOrg.chr_opts
	o.generation = 0
	o.dead = false
	return o
end

local function kill(o)
	o.dead = true
	o.homologs = nil
	return o
end

local function checkpoint(org)
	if org.dead or not org.homologs then
		return { repaired = 0, dead = 0, structural = 0, anomaly = 0 }
	end
	local rep = { repaired = 0, dead = 0, structural = 0, anomaly = 0 }
	for i = 1, #org.chr_ids do
		local a3, b2, r = cell.checkpoint_pair(org.homologs[1][i], org.homologs[2][i])
		org.homologs[1][i], org.homologs[2][i] = a3, b2
		rep.repaired = rep.repaired + r.repaired
		rep.dead = rep.dead + r.dead
		rep.structural = rep.structural + r.structural
		rep.anomaly = rep.anomaly + r.anomaly
	end
	return rep
end

local function read_checked(org)
	local rep = checkpoint(org)
	local data = {}
	for i = 1, #org.chr_ids do
		local d
		for h = 1, 2 do
			local rec = chromosome.parse(org.homologs[h][i])
			if rec and rec.cen_ok then
				d = chromosome.read(rec)
				if d then break end
			end
		end
		if not d then return nil, org.chr_ids[i], rep end
		data[org.chr_ids[i]] = d
	end
	return data, nil, rep
end

local function read(org)
	if org.dead or not org.homologs then
		if org.stem_source and renew(org, org.stem_source) then
			local data = read_checked(org)
			if data then return data, { renewed = true } end
		end
		return nil, "organism dead"
	end
	local data, broken_chr, rep = read_checked(org)
	if data then return data, rep end
	if org.stem_source and renew(org, org.stem_source) then
		local data2 = read_checked(org)
		if data2 then
			rep.renewed = true
			return data2, rep
		end
	end
	return nil, "chromosome " .. tostring(broken_chr) .. " dead", rep
end

local function replicate(org)
	if org.dead or not org.homologs then return nil, "organism dead" end
	checkpoint(org)
	if org.generation + 1 > org.max_gen then return nil, "senescent" end
	org.generation = org.generation + 1
	for i = 1, #org.chr_ids do
		local s = chromosome.set_generation(org.homologs[1][i], org.generation)
			or chromosome.set_generation(org.homologs[2][i], org.generation)
		if not s then return nil, "centromere damaged" end
		org.homologs[1][i], org.homologs[2][i] = s, s
	end
	return org
end

local function mitosis(org)
	if org.dead or not org.homologs then return nil, "organism dead" end
	checkpoint(org)
	if org.generation + 1 > org.max_gen then return nil, "senescent" end
	org.generation = org.generation + 1
	local gA, gB = {}, {}
	for i = 1, #org.chr_ids do
		local s = chromosome.set_generation(org.homologs[1][i], org.generation)
			or chromosome.set_generation(org.homologs[2][i], org.generation)
		if not s then return nil, "centromere damaged" end
		org.homologs[1][i], org.homologs[2][i] = s, s
		gA[i], gB[i] = s, s
	end
	return {
		chr_ids = org.chr_ids,
		chr_opts = org.chr_opts,
		homologs = { gA, gB },
		generation = org.generation,
		max_gen = org.max_gen,
		stem = false,
		dead = false,
		stem_source = org.stem_source,
	}
end

local function damage(org, count, seed)
	if org.dead or not org.homologs then return nil, "organism dead" end
	for i = 1, #org.chr_ids do
		local h1 = cell.damage_strand(org.homologs[1][i], count, (seed or 1) + i)
		local h2 = cell.damage_strand(org.homologs[2][i], count, (seed or 1) + 1000 + i)
		org.homologs[1][i], org.homologs[2][i] = h1, h2
	end
	return org
end

local function mutate(org, count, seed)
	if org.dead or not org.homologs then return nil, "organism dead" end
	local rng = cell.rng(seed)
	local ci = 1 + rng() % #org.chr_ids
	local cid = org.chr_ids[ci]
	local data
	for h = 1, 2 do
		local rec = chromosome.parse(org.homologs[h][ci])
		if rec and rec.cen_ok then
			data = chromosome.read(rec)
			if data then break end
		end
	end
	if not data then return nil, "chromosome " .. cid .. " unreadable" end
	local t = { byte(data, 1, #data) }
	if #t == 0 or count == nil or count <= 0 then
		return { chr = cid, data = data }
	end
	for _ = 1, count do
		local bi = 1 + rng() % #t
		t[bi] = bxor(t[bi], lshift(1, rng() % 8))
	end
	local mutated = bytes_to_string(t)
	local strand, err = chromosome.encode(cid, mutated, org.chr_opts[ci])
	if not strand then return nil, "re-encode failed: " .. tostring(err) end
	org.homologs[1][ci], org.homologs[2][ci] = strand, strand
	return { chr = cid, data = mutated }
end

local function cross(pa, pb, seed)
	if type(pa) ~= "table" or type(pb) ~= "table"
		or not pa.homologs or not pb.homologs then
		return nil, "expected organisms"
	end
	if #pa.chr_ids ~= #pb.chr_ids then return nil, "incompatible genomes" end
	local rng = cell.rng(seed)
	local ids, opts_t = {}, {}
	local gA, gB = {}, {}
	for i = 1, #pa.chr_ids do
		if pa.chr_ids[i] ~= pb.chr_ids[i] then return nil, "incompatible genomes" end
		if not opts_equal(pa.chr_opts[i], pb.chr_opts[i]) then
			return nil, "incompatible chromosome opts"
		end
		local ga, gb, ngenesA, ngenesB = {}, {}, nil, nil
		for h = 1, 2 do
			local recA = chromosome.parse(pa.homologs[h][i])
			if recA and recA.cen_ok then
				ngenesA = recA.ngenes
				for _, g in ipairs(recA.genes) do
					if g.crc_ok and ga[g.id] == nil then ga[g.id] = g.data end
				end
			end
			local recB = chromosome.parse(pb.homologs[h][i])
			if recB and recB.cen_ok then
				ngenesB = recB.ngenes
				for _, g in ipairs(recB.genes) do
					if g.crc_ok and gb[g.id] == nil then gb[g.id] = g.data end
				end
			end
		end
		if not ngenesA or not ngenesB then return nil, "parent chromosome broken" end
		if ngenesA ~= ngenesB then return nil, "incompatible gene counts on chr " .. pa.chr_ids[i] end

		-- Build one balanced recombinant mosaic chromosome for both homologs
		local parts = {}
		for gid = 0, ngenesA - 1 do
			local useB = rng() % 2 == 1
			local d = useB and (gb[gid] or ga[gid]) or (ga[gid] or gb[gid])
			if not d then return nil, "parent gene " .. gid .. " missing" end
			parts[#parts + 1] = d
		end
		local strand, err = chromosome.encode(pa.chr_ids[i], concat(parts), pa.chr_opts[i])
		if not strand then return nil, "cross re-encode failed: " .. tostring(err) end
		gA[i], gB[i] = strand, strand
		ids[i] = pa.chr_ids[i]
		opts_t[i] = pa.chr_opts[i]
	end
	return {
		chr_ids = ids,
		chr_opts = opts_t,
		homologs = { gA, gB },
		generation = 0,
		max_gen = pa.max_gen,
		stem = false,
		dead = false,
		stem_source = nil,
	}
end

local function serialize(org)
	if org.dead or not org.homologs then return nil, "organism dead" end
	local t = { MAGIC,
		char(band(rshift(org.generation, 8), 255), band(org.generation, 255)),
		char(#org.chr_ids) }
	for i = 1, #org.chr_ids do
		local strand = org.homologs[1][i]
		local o = org.chr_opts[i] or {}
		local mode = (o.mode or "dense") == "codon" and 1 or 0
		local h = (o.h or 3) - 1
		local gene_raw = o.gene_raw or 1024
		t[#t + 1] = char(org.chr_ids[i],
			band(rshift(gene_raw, 8), 255), band(gene_raw, 255),
			band(bor(mode, lshift(h, 1)), 255),
			band(o.units or 4, 255), band(o.flags or 0, 255))
		t[#t + 1] = char(band(rshift(#strand, 24), 255), band(rshift(#strand, 16), 255),
			band(rshift(#strand, 8), 255), band(#strand, 255))
		t[#t + 1] = strand
	end
	return concat(t)
end

-- VIV14 / VIV14N: versioned container with both homologs and max_gen.
-- VIV14N additionally stores the parity count per chromosome; when a
-- VIV14 container carries CEN2 strands the count is inferred from them.
local function serialize_rev(org, magic, with_parity, with_primers)
	if org.dead or not org.homologs then return nil, "organism dead" end
	if org.generation < 0 or org.generation > 65535 then
		return nil, "generation out of range"
	end
	if org.max_gen < 0 or org.max_gen > 65535 then
		return nil, "max_gen out of range"
	end
	if with_primers then
		for i = 1, #org.chr_ids do
			local p = (org.chr_opts[i] or {}).primer or 0
			if type(p) ~= "number" or p ~= floor(p) or p < 0 or p > 65535 then
				return nil, "primer out of range"
			end
		end
	end
	local t = { magic, char(1),  -- flags: diploid pair stored
		char(band(rshift(org.generation, 8), 255), band(org.generation, 255)),
		char(#org.chr_ids),
		char(band(rshift(org.max_gen, 8), 255), band(org.max_gen, 255)) }
	for i = 1, #org.chr_ids do
		local s0, s1 = org.homologs[1][i], org.homologs[2][i]
		local o = org.chr_opts[i] or {}
		local mode = (o.mode or "dense") == "codon" and 1 or 0
		local h = (o.h or 3) - 1
		local gene_raw = o.gene_raw or 1024
		local head = char(org.chr_ids[i],
			band(rshift(gene_raw, 8), 255), band(gene_raw, 255),
			band(bor(mode, lshift(h, 1)), 255),
			band(o.units or 4, 255), band(o.flags or 0, 255))
		if with_parity then
			head = head .. char(band(o.parity or 0, 255))
		end
		if with_primers then
			head = head .. char(band(rshift(o.primer or 0, 8), 255),
				band(o.primer or 0, 255))
		end
		t[#t + 1] = head
		t[#t + 1] = char(band(rshift(#s0, 24), 255), band(rshift(#s0, 16), 255),
			band(rshift(#s0, 8), 255), band(#s0, 255))
		t[#t + 1] = s0
		t[#t + 1] = char(band(rshift(#s1, 24), 255), band(rshift(#s1, 16), 255),
			band(rshift(#s1, 8), 255), band(#s1, 255))
		t[#t + 1] = s1
	end
	return concat(t)
end

local function serialize14(org)
	return serialize_rev(org, MAGIC14, false, false)
end

local function serialize14n(org)
	return serialize_rev(org, MAGIC14N, true, false)
end

local function serialize14nb(org)
	return serialize_rev(org, MAGIC_BANSHEE, true, true)
end

local function container_version(s)
	if type(s) ~= "string" then return 0 end
	-- "VIV14NB4NSH33" is reserved (not implemented yet)
	if #s >= 13 and sub(s, 1, 13) == MAGIC_BANSHEE then return 143 end
	if #s >= 6 and sub(s, 1, 6) == MAGIC14N then return 141 end
	if #s >= 5 and sub(s, 1, 5) == MAGIC14 then return 14 end
	if #s >= 4 and sub(s, 1, 4) == MAGIC then return 1 end
	return 0
end

local function deserialize_viv1(s)
	if type(s) ~= "string" or #s < 7 then return nil, "bad format" end
	if sub(s, 1, 4) ~= MAGIC then return nil, "bad magic" end
	local generation = byte(s, 5) * 256 + byte(s, 6)
	if generation > 65535 then return nil, "bad generation" end
	local nchr = byte(s, 7)
	if nchr < 1 then return nil, "bad chromosome count" end
	local pos = 8
	local ids, opts_t, gA, gB = {}, {}, {}, {}
	for i = 1, nchr do
		if pos + 9 > #s then return nil, "truncated" end
		local cid = byte(s, pos)
		local gene_raw = byte(s, pos + 1) * 256 + byte(s, pos + 2)
		local mode_h = byte(s, pos + 3)
		local mode = (band(mode_h, 1) == 1) and "codon" or "dense"
		local h = band(rshift(mode_h, 1), 15) + 1
		local units = byte(s, pos + 4)
		local flags = byte(s, pos + 5)
		local slen = byte(s, pos + 6) * 0x1000000 + byte(s, pos + 7) * 0x10000
			+ byte(s, pos + 8) * 0x100 + byte(s, pos + 9)
		pos = pos + 10
		if pos + slen - 1 > #s then return nil, "truncated" end
		local strand = sub(s, pos, pos + slen - 1)
		pos = pos + slen
		if gene_raw < 16 or h < 3 or h > 12 or units < 1 then
			return nil, "corrupt chromosome options"
		end
		local rec = chromosome.parse(strand)
		if not rec or not rec.cen_ok or rec.id ~= cid then
			return nil, "corrupt chromosome " .. i
		end
		ids[i] = cid
		opts_t[i] = { gene_raw = gene_raw, mode = mode, h = h, units = units, flags = flags }
		gA[i], gB[i] = strand, strand
	end
	return {
		chr_ids = ids,
		chr_opts = opts_t,
		homologs = { gA, gB },
		generation = generation,
		max_gen = 60,
		stem = false,
		dead = false,
		stem_source = nil,
	}
end

local function be32(s, p)
	return byte(s, p) * 0x1000000 + byte(s, p + 1) * 0x10000
		+ byte(s, p + 2) * 0x100 + byte(s, p + 3)
end

local function deserialize_rev(s, with_parity, with_primers)
	local mlen = with_primers and 13 or (with_parity and 6 or 5)
	local hdr = mlen + 6
	local chdr = 6 + (with_parity and 1 or 0) + (with_primers and 2 or 0)
	if #s < hdr then return nil, "bad format" end
	if sub(s, 1, mlen) ~= (with_primers and MAGIC_BANSHEE
		or (with_parity and MAGIC14N or MAGIC14)) then
		return nil, "bad magic"
	end
	if byte(s, mlen + 1) ~= 1 then return nil, "unsupported flags" end
	local generation = byte(s, mlen + 2) * 256 + byte(s, mlen + 3)
	local nchr = byte(s, mlen + 4)
	local max_gen = byte(s, mlen + 5) * 256 + byte(s, mlen + 6)
	if nchr < 1 then return nil, "bad chromosome count" end
	local pos = hdr + 1
	local ids, opts_t, gA, gB = {}, {}, {}, {}
	for i = 1, nchr do
		if pos + chdr + 4 - 1 > #s then return nil, "truncated" end
		local cid = byte(s, pos)
		local gene_raw = byte(s, pos + 1) * 256 + byte(s, pos + 2)
		local mode_h = byte(s, pos + 3)
		local mode = (band(mode_h, 1) == 1) and "codon" or "dense"
		local h = band(rshift(mode_h, 1), 15) + 1
		local units = byte(s, pos + 4)
		local flags = byte(s, pos + 5)
		local parity_byte = with_parity and byte(s, pos + 6) or 0
		local primer_byte = with_primers and (byte(s, pos + 7) * 256 + byte(s, pos + 8)) or 0
		pos = pos + chdr
		if gene_raw < 16 or h < 3 or h > 12 or units < 1 then
			return nil, "corrupt chromosome options"
		end
		local sl0 = be32(s, pos)
		pos = pos + 4
		if pos + sl0 - 1 > #s then return nil, "truncated" end
		local s0 = sub(s, pos, pos + sl0 - 1)
		pos = pos + sl0
		local sl1 = be32(s, pos)
		pos = pos + 4
		if pos + sl1 - 1 > #s then return nil, "truncated" end
		local s1 = sub(s, pos, pos + sl1 - 1)
		pos = pos + sl1
		-- both homologs must parse; at least one centromere must be intact
		local r0 = chromosome.parse(s0)
		local r1 = chromosome.parse(s1)
		if not r0 or not r1 then return nil, "corrupt chromosome" end
		local cen0 = r0.cen_ok and r0.id == cid
		local cen1 = r1.cen_ok and r1.id == cid
		if not cen0 and not cen1 then return nil, "corrupt chromosome" end
		ids[i] = cid
		opts_t[i] = {
			gene_raw = gene_raw, mode = mode, h = h, units = units, flags = flags,
			-- VIV14N stores the parity count; VIV14 infers it from the strand
			parity = with_parity and parity_byte or r0.parity,
			-- Banshee stores the barcode; earlier containers infer it
			primer = with_primers and primer_byte or r0.primer,
		}
		gA[i], gB[i] = s0, s1
	end
	return {
		chr_ids = ids,
		chr_opts = opts_t,
		homologs = { gA, gB },
		generation = generation,
		max_gen = max_gen ~= 0 and max_gen or 60,
		stem = false,
		dead = false,
		stem_source = nil,
	}
end

local function deserialize(s)
	if type(s) ~= "string" then return nil, "bad format" end
	local rev = container_version(s)
	if rev == 143 then return deserialize_rev(s, true, true) end
	if rev == 141 then return deserialize_rev(s, true, false) end
	if rev == 14 then
		local o = deserialize_rev(s, false, false)
		if o then return o end
		-- VIV1's generation high byte can be '4', so a VIV1 file may look
		-- like a VIV14 one: fall back when the VIV14 structure is invalid
	end
	return deserialize_viv1(s)
end

local function maintain(organisms, stemOrg)
	local report = { checked = #organisms, repaired = 0, renewed = 0, dead = 0 }
	for i = 1, #organisms do
		local o = organisms[i]
		if o.dead or not o.homologs then
			if stemOrg and renew(o, stemOrg) then
				report.renewed = report.renewed + 1
			else
				kill(o)
				report.dead = report.dead + 1
			end
		else
			local data, _, rep = read_checked(o)
			report.repaired = report.repaired + rep.repaired + rep.structural
			if not data then
				if stemOrg and renew(o, stemOrg) then
					report.renewed = report.renewed + 1
				else
					kill(o)
					report.dead = report.dead + 1
				end
			elseif o.generation >= o.max_gen and stemOrg then
				renew(o, stemOrg)
				report.renewed = report.renewed + 1
			end
		end
	end
	return report
end

return {
	new = new,
	stem = stem,
	attach_stem = attach_stem,
	read = read,
	checkpoint = checkpoint,
	replicate = replicate,
	mitosis = mitosis,
	renew = renew,
	kill = kill,
	damage = damage,
	mutate = mutate,
	cross = cross,
	serialize = serialize,
	serialize14 = serialize14,
	serialize14n = serialize14n,
	serialize14nb = serialize14nb,
	container_version = container_version,
	deserialize = deserialize,
	maintain = maintain,
}
