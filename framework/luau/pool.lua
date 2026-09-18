--!native
local channel = require("./channel")
local genome = require("./genome")
local sub = string.sub

local TWO32 = 4294967296
local CROSS_MUL = 0x9E3779B9

local function new()
	return { ids = {}, strands = {} }
end

local function add_gene(p, id, strand)
	p.ids[#p.ids + 1] = id
	p.strands[#p.strands + 1] = strand
	return true
end

local function add_chromosome(p, strand)
	local genes, err = genome.scan(strand)
	if not genes then return nil, err end
	for i = 1, #genes do
		local g = genes[i]
		if g.crc_ok then
			add_gene(p, g.id, sub(strand, g.offset, g.offset + g.size - 1))
		end
	end
	return true
end

local function amplify(p, id, opts)
	opts = opts or {}
	local p_access = opts.p_access or 0.0
	local p_cross = opts.p_cross or 0.0
	local seed = opts.seed or 0
	local ch = opts.ch or {}
	local p_primer = opts.p_primer or 0.0
	local coverage = opts.coverage or 1
	local n = #p.ids
	if n == 0 then return nil, "empty pool" end
	if not (p_access >= 0.0 and p_access <= 1.0)
		or not (p_cross >= 0.0 and p_cross <= 1.0)
		or not (p_primer >= 0.0 and p_primer <= 1.0) then
		return nil, "invalid probability"
	end
	local target
	for i = 1, n do
		if p.ids[i] == id then
			target = i
			break
		end
	end
	if not target then return nil, "target not in pool" end
	local rng = channel.prng_new(seed)
	if channel.prng_chance(rng, p_access) then
		return { read = { strand = "", bases = 0, dropped = 1 }, id = -1 }
	end
	if p_primer > 0.0 and channel.prng_chance(rng, p_primer) then
		return { read = { strand = "", bases = 0, dropped = 1 }, id = -1 }
	end
	local pick = target
	if n > 1 and channel.prng_chance(rng, p_cross) then
		pick = channel.prng_next(rng) % n + 1
		if p.ids[pick] == id then pick = pick % n + 1 end
	end
	local c = {}
	for k, v in pairs(ch) do c[k] = v end
	c.seed = (seed + (pick - 1) * CROSS_MUL % TWO32 + 1) % TWO32
	local rd, err = channel.consensus_read(p.strands[pick], { ch = c, coverage = coverage })
	if not rd then return nil, err end
	return { read = rd, id = p.ids[pick] }
end

return {
	new = new,
	add_gene = add_gene,
	add_chromosome = add_chromosome,
	amplify = amplify,
}
