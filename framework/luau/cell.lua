--!native
-- cell.lua  (phases 3-4 -- the diploid cell and self-renewal)
--
-- Phase 3 (diploidy + repair): every chromosome is kept as a HOMOLOG
-- PAIR (maternal/paternal copies), like a real diploid genome.
-- checkpoint() compares the homologs and repairs:
--   * a damaged gene is excised and resynthesized from the healthy
--     homolog at the same locus (homologous recombination -- damage
--     flips bytes but never shifts the layout, so loci align)
--   * structural damage (telomeres / centromere) is spliced from the
--     healthy homolog; a destroyed centromere triggers a full rebuild
--   * a locus damaged in BOTH homologs is unrecoverable internally
--     (dead gene) -- phase 4 replaces the cell from a stem cell
--
-- Phase 4 (the cell): replicate with proofreading, mitosis, the
-- Hayflick limit (senescence), apoptosis, and self-renewal: a stem
-- cell serves as the readonly germ-line reference from which damaged
-- or aged cells are regenerated (generation resets to 0).
--
-- Requires bit32 (Luau, Lua 5.2+), chromosome, genome, dna.

local chromosome = require("./chromosome")

local byte, char, concat, sub =
	string.byte, string.char, table.concat, string.sub
local floor = math.floor

local band, bxor, lshift, rshift, bor =
	bit32.band, bit32.bxor, bit32.lshift, bit32.rshift, bit32.bor

-- deterministic xorshift32 rng for reproducible damage
local function make_rng(seed)
	local state = (seed or 1) % 4294967296
	if state == 0 then state = 1 end
	return function()
		state = bxor(state, lshift(state, 13))
		state = bxor(state, rshift(state, 17))
		state = bxor(state, lshift(state, 5))
		return state
	end
end

local function damage_strand(strand, count, seed)
	local s = (type(strand) == "string") and strand or nil
	if s == nil then return nil, "expected string" end
	local n = #s * 4
	if n < 1 then return "" end
	if count == nil or count <= 0 then return s .. "" end  -- unchanged strand
	local rng = make_rng(seed)
	local edits = {}
	for _ = 1, count do
		local idx = rng() % n
		local bytei = floor(idx / 4) + 1
		local sh = 2 * (3 - idx % 4)
		local d = band(rshift(byte(s, bytei), sh), 3)
		edits[idx] = (d + 1 + rng() % 3) % 4
	end
	local t = {}
	for i = 1, #s do
		local b = byte(s, i)
		for j = 0, 3 do
			local idx = (i - 1) * 4 + j
			if edits[idx] then
				local sh = 2 * (3 - j)
				b = bor(band(b, bxor(255, lshift(3, sh))), lshift(edits[idx], sh))
			end
		end
		t[i] = char(b)
	end
	return concat(t)
end

-- phase 3: repair `dst` using homolog `src` as the template
local function repair_homolog(dst, src)
	local report = { repaired = 0, dead = 0, structural = 0, anomaly = 0 }
	local recS = chromosome.parse(src)
	if not recS or not recS.cen_ok then
		report.dead = report.dead + 1
		return dst, report
	end
	local tb = recS.telomere_bytes
	local recD = chromosome.parse(dst)
	if not recD or not recD.cen_ok then
		-- centromere destroyed in dst: the head region (left telomere +
		-- primer site + centromere) lives at fixed positions -- splice it
		-- wholesale, with the source layout
		local cen_bytes = (recS.cen_version == 2) and chromosome.CEN2BYTES
			or chromosome.CENBYTES
		local head = tb + (recS.primer_bytes or 0) + cen_bytes
		dst = sub(src, 1, head) .. sub(dst, head + 1)
		report.structural = report.structural + 1
		recD = chromosome.parse(dst)
	end
	if recD and not recD.telomere_ok and recS.telomere_ok then
		local fixed = false
		if sub(dst, 1, tb) ~= sub(src, 1, tb) then
			dst = sub(src, 1, tb) .. sub(dst, tb + 1)
			fixed = true
		end
		if sub(dst, #dst - tb + 1) ~= sub(src, #src - tb + 1) then
			dst = sub(dst, 1, #dst - tb) .. sub(src, #src - tb + 1)
			fixed = true
		end
		if fixed then
			report.structural = report.structural + 1
			recD = chromosome.parse(dst)
		end
	end
	if not recD or not recD.cen_ok then
		report.dead = report.dead + 1
		return dst, report
	end
	local dmap = {}
	for i = 1, #recD.genes do
		local g = recD.genes[i]
		if not dmap[g.id] then dmap[g.id] = g end
	end
	for i = 1, #recS.genes do
		local g = recS.genes[i]
		local d = dmap[g.id]
		local need = false
		if not g.crc_ok then
			-- dead only if broken in BOTH homologs: a healthy dst gene
			-- will repair src on the next checkpoint pass
			if not d or not d.crc_ok then
				report.dead = report.dead + 1
			end
		elseif not d then
			need = true                    -- gene lost in dst (promoter hit)
		elseif not d.crc_ok or d.data ~= g.data then
			need = true                    -- damaged or altered in dst
		end
		if need and g.crc_ok then
			-- the true locus always equals the template's offset: point
			-- damage never shifts the layout, and a phantom promoter in
			-- dst must never redirect the splice
			dst = sub(dst, 1, g.offset - 1)
				.. sub(src, g.offset, g.offset + g.size - 1)
				.. sub(dst, g.offset + g.size)
			report.repaired = report.repaired + 1
			dmap[g.id] = { id = g.id, offset = g.offset, size = g.size,
				crc_ok = true, data = g.data }
		end
	end
	if #recD.genes > #recS.genes then
		report.anomaly = report.anomaly + (#recD.genes - #recS.genes)
	end
	return dst, report
end

-- phase 3: diploid checkpoint -- symmetric cascade of repairs for one
-- chromosome pair: each homolog repairs the other, and the first one
-- gets a third pass so that a homolog healed on the second pass can
-- rescue the first one (e.g. broken gene here + broken centromere
-- there). "dead" counts both-broken encounters during the passes;
-- actual recoverability is judged by the readability checks afterwards
local function checkpoint_pair(h1, h2)
	local a2, rep1 = repair_homolog(h1, h2)
	local b2, rep2 = repair_homolog(h2, a2)
	local a3, rep3 = repair_homolog(a2, b2)
	return a3, b2, {
		repaired = rep1.repaired + rep2.repaired + rep3.repaired,
		dead = rep1.dead + rep2.dead + rep3.dead,
		structural = rep1.structural + rep2.structural + rep3.structural,
		anomaly = rep1.anomaly + rep2.anomaly + rep3.anomaly,
	}
end

local function checkpoint(c)
	if c.dead or not c.homologs then
		return { repaired = 0, dead = 0, structural = 0, anomaly = 0 }
	end
	local a3, b2, rep = checkpoint_pair(c.homologs[1], c.homologs[2])
	c.homologs[1], c.homologs[2] = a3, b2
	return rep
end

-- phase 4: the cell
local function new(chr_id, data, opts)
	opts = opts or {}
	local strand, err = chromosome.encode(chr_id, data, opts)
	if not strand then return nil, err or "chromosome encode failed" end
	return {
		chr_id = chr_id,
		homologs = { strand, strand },
		generation = 0,
		max_gen = opts.max_gen or 60,
		stem = false,
		dead = false,
		stem_source = nil,
	}
end

local function stem(chr_id, data, opts)
	local c, err = new(chr_id, data, opts)
	if not c then return nil, err end
	c.stem = true
	return c
end

local function attach_stem(c, stemc)
	c.stem_source = stemc
	return c
end

local function renew(c, stemc)
	if not stemc or not stemc.homologs then return nil, "no stem cell" end
	local rec = chromosome.parse(stemc.homologs[1])
	if not rec or not rec.cen_ok then return nil, "stem cell damaged" end
	c.homologs = { stemc.homologs[1], stemc.homologs[2] }
	c.generation = 0
	c.dead = false
	return c
end

local function kill(c)
	c.dead = true
	c.homologs = nil
	return c
end

local function damage(c, count, seed)
	if c.dead or not c.homologs then return nil, "cell is dead" end
	local a = damage_strand(c.homologs[1], count, seed)
	local b = damage_strand(c.homologs[2], count, (seed or 1) + 1)
	c.homologs = { a, b }
	return c
end

-- reads through the repair machinery; falls back to the homolog, then
-- to self-renewal from the stem cell
local function read(c)
	if c.dead or not c.homologs then
		if c.stem_source and renew(c, c.stem_source) then
			local rec = chromosome.parse(c.homologs[1])
			local data = rec and chromosome.read(rec)
			if data then return data, { renewed = true } end
		end
		return nil, "cell is dead"
	end
	local rep = checkpoint(c)
	for i = 1, 2 do
		local rec = chromosome.parse(c.homologs[i])
		if rec and rec.cen_ok then
			local data = chromosome.read(rec)
			if data then return data, rep end
		end
	end
	if c.stem_source and renew(c, c.stem_source) then
		local rec = chromosome.parse(c.homologs[1])
		local data = rec and chromosome.read(rec)
		if data then
			rep.renewed = true
			return data, rep
		end
	end
	return nil, "cell is dead", rep
end

-- proofread replication: checkpoint, then fresh homologs with the
-- generation counter advanced
local function replicate(c)
	if c.dead or not c.homologs then return nil, "cell is dead" end
	checkpoint(c)
	if c.generation + 1 > c.max_gen then return nil, "senescent" end
	c.generation = c.generation + 1
	local s = chromosome.set_generation(c.homologs[1], c.generation)
		or chromosome.set_generation(c.homologs[2], c.generation)
	if not s then return nil, "centromere damaged" end
	c.homologs = { s, s }
	return c
end

-- mitosis: the mother cell also ages (a division consumes the mother):
-- both the parent and the identical daughter advance by one generation
local function mitosis(c)
	if c.dead or not c.homologs then return nil, "cell is dead" end
	checkpoint(c)
	if c.generation + 1 > c.max_gen then return nil, "senescent" end
	c.generation = c.generation + 1
	local s = chromosome.set_generation(c.homologs[1], c.generation)
		or chromosome.set_generation(c.homologs[2], c.generation)
	if not s then return nil, "centromere damaged" end
	c.homologs = { s, s }
	return {
		chr_id = c.chr_id,
		homologs = { s, s },
		generation = c.generation,
		max_gen = c.max_gen,
		stem = false,
		dead = false,
		stem_source = c.stem_source,
	}
end

-- the self-renewal loop: checkpoint the population, replace dead,
-- damaged-beyond-repair or aged cells from the stem cell
local function maintain(cells, stemc)
	local report = { checked = #cells, repaired = 0, renewed = 0, dead = 0 }
	for i = 1, #cells do
		local c = cells[i]
		if c.dead or not c.homologs then
			if stemc and renew(c, stemc) then
				report.renewed = report.renewed + 1
			else
				report.dead = report.dead + 1
			end
		else
			local rep = checkpoint(c)
			report.repaired = report.repaired + rep.repaired + rep.structural
			local data
			for h = 1, 2 do
				local rec = chromosome.parse(c.homologs[h])
				if rec and rec.cen_ok then
					data = chromosome.read(rec)
					if data then break end
				end
			end
			if not data then
				if stemc and renew(c, stemc) then
					report.renewed = report.renewed + 1
				else
					kill(c)  -- mark dead so the next cycle skips it
					report.dead = report.dead + 1
				end
			elseif c.generation >= c.max_gen and stemc then
				renew(c, stemc)  -- Hayflick limit: regenerate from the niche
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
	checkpoint_pair = checkpoint_pair,
	repair_homolog = repair_homolog,
	rng = make_rng,
	replicate = replicate,
	mitosis = mitosis,
	renew = renew,
	kill = kill,
	damage = damage,
	damage_strand = damage_strand,
	maintain = maintain,
}
