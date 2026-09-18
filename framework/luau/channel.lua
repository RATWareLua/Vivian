--!native
local band, bxor, rshift, lshift, bor =
	bit32.band, bit32.bxor, bit32.rshift, bit32.lshift, bit32.bor
local floor = math.floor
local byte, char, concat = string.byte, string.char, table.concat

local TWO32 = 4294967296
local TWO16 = 65536
local ADD_HI, ADD_LO = 0x9E3779B9, 0x7F4A7C15
local M1_HI, M1_LO = 0xBF58476D, 0x1CE4E5B9
local M2_HI, M2_LO = 0x94D049BB, 0x133111EB

local function add64(a, b)
	local lo = a.lo + b.lo
	local hi = a.hi + b.hi
	if lo >= TWO32 then
		lo = lo - TWO32
		hi = hi + 1
	end
	return { hi = hi % TWO32, lo = lo }
end

local function xor64(a, b)
	return { hi = bxor(a.hi, b.hi), lo = bxor(a.lo, b.lo) }
end

local function shr64(a, n)
	if n == 0 then
		return { hi = a.hi, lo = a.lo }
	elseif n < 32 then
		return { hi = rshift(a.hi, n), lo = bor(rshift(a.lo, n), lshift(a.hi, 32 - n)) }
	elseif n == 32 then
		return { hi = 0, lo = a.hi }
	end
	return { hi = 0, lo = rshift(a.hi, n - 32) }
end

local function mul64(a, mhi, mlo)
	local a0, a1 = band(a.lo, 0xFFFF), rshift(a.lo, 16)
	local a2, a3 = band(a.hi, 0xFFFF), rshift(a.hi, 16)
	local b0, b1 = band(mlo, 0xFFFF), rshift(mlo, 16)
	local b2, b3 = band(mhi, 0xFFFF), rshift(mhi, 16)
	local p0 = a0 * b0
	local p1 = a0 * b1 + a1 * b0
	local p2 = a0 * b2 + a1 * b1 + a2 * b0
	local p3 = a0 * b3 + a1 * b2 + a2 * b1 + a3 * b0
	local r0 = p0 % TWO16
	local c = floor(p0 / TWO16)
	local t1 = p1 + c
	local r1 = t1 % TWO16
	c = floor(t1 / TWO16)
	local t2 = p2 + c
	local r2 = t2 % TWO16
	c = floor(t2 / TWO16)
	local t3 = p3 + c
	local r3 = t3 % TWO16
	return { hi = r2 + r3 * TWO16, lo = r0 + r1 * TWO16 }
end

local function prng_new(seed)
	return { state = { hi = 0, lo = seed % TWO32 } }
end

local function prng_mix(p)
	p.state = add64(p.state, { hi = ADD_HI, lo = ADD_LO })
	local z = p.state
	z = xor64(z, shr64(z, 30))
	z = mul64(z, M1_HI, M1_LO)
	z = xor64(z, shr64(z, 27))
	z = mul64(z, M2_HI, M2_LO)
	return xor64(z, shr64(z, 31))
end

local function prng_next(p)
	return prng_mix(p).hi
end

local function prng_next64(p)
	return prng_mix(p)
end

local function prng_chance(p, prob)
	if prob <= 0.0 then return false end
	if prob >= 1.0 then return true end
	local thr = floor(prob * 18446744073709551616.0)
	local r = prng_mix(p)
	local thi = floor(thr / TWO32)
	if r.hi ~= thi then return r.hi < thi end
	return r.lo < thr % TWO32
end

local function digit_at(s, i)
	return band(rshift(byte(s, floor(i / 4) + 1), 2 * (3 - i % 4)), 3)
end

local function pack_digits(digits, m)
	if m == 0 then return "" end
	local bytes = floor((m + 3) / 4)
	local out = {}
	for i = 1, bytes do out[i] = 0 end
	for i = 0, m - 1 do
		local k = floor(i / 4) + 1
		out[k] = bor(out[k], lshift(digits[i + 1], 2 * (3 - i % 4)))
	end
	local t = {}
	for i = 1, bytes do t[i] = char(out[i]) end
	return concat(t)
end

local function channel_read_impl(strand, opts, soft)
	if type(strand) ~= "string" then return nil, "expected string" end
	opts = opts or {}
	local p_sub = opts.p_sub or 0.0
	local p_ins = opts.p_ins or 0.0
	local p_del = opts.p_del or 0.0
	local p_drop = opts.p_drop or 0.0
	local seed = opts.seed or 0
	local p_sub_gc = opts.p_sub_gc or 0.0
	local p_sub_hp = opts.p_sub_hp or 0.0
	local p_trunc = opts.p_trunc or 0.0
	local p_burst = opts.p_burst or 0.0
	local burst_len = opts.burst_len or 0
	local p_burst_del = opts.p_burst_del or 0.0
	local ps = { p_sub, p_ins, p_del, p_drop, p_sub_gc, p_sub_hp, p_trunc, p_burst, p_burst_del }
	for i = 1, 9 do
		if not (ps[i] >= 0.0 and ps[i] <= 1.0) then
			return nil, "invalid probability"
		end
	end
	local rng = prng_new(seed)
	if prng_chance(rng, p_drop) then
		return { strand = "", bases = 0, dropped = 1 }
	end
	local n = #strand * 4
	local limit = 0
	if p_trunc > 0.0 and n >= 2 and prng_chance(rng, p_trunc) then
		limit = 1 + prng_next(rng) % (n - 1)
	end
	local digits = {}
	local qual = soft and {} or nil
	local m = 0
	local burst_rem = 0
	for i = 0, n - 1 do
		if limit ~= 0 and m >= limit then break end
		local d = digit_at(strand, i)
		local bursted = false
		if p_burst > 0.0 then
			if burst_rem > 0 then
				bursted = true
				burst_rem = burst_rem - 1
			elseif prng_chance(rng, p_burst) then
				burst_rem = (burst_len ~= 0) and burst_len or 8
				bursted = true
			end
		end
		if not prng_chance(rng, p_del) then
			local psub = p_sub
			local killed = false
			if bursted then
				if p_burst_del > 0.0 and prng_chance(rng, p_burst_del) then
					killed = true
				else
					psub = 1.0
				end
			else
				if p_sub_gc > 0.0 and (d == 1 or d == 2) then
					psub = 1.0 - (1.0 - psub) * (1.0 - p_sub_gc)
				end
				if p_sub_hp > 0.0 and i > 0 and d == digit_at(strand, i - 1) then
					psub = 1.0 - (1.0 - psub) * (1.0 - p_sub_hp)
				end
			end
			if not killed then
				local sub = prng_chance(rng, psub)
				if sub then
					local k = prng_next(rng) % 3
					d = (d + 1 + k) % 4
				end
				if limit == 0 or m < limit then
					m = m + 1
					digits[m] = d
					if qual then qual[m] = sub and 4 or 60 end
				end
			end
		end
		if prng_chance(rng, p_ins) then
			local ins = prng_next(rng) % 4
			if limit == 0 or m < limit then
				m = m + 1
				digits[m] = ins
				if qual then qual[m] = 4 end
			end
		end
	end
	if soft then
		return { strand = pack_digits(digits, m), bases = m, dropped = 0, qual = qual }
	end
	return { strand = pack_digits(digits, m), bases = m, dropped = 0 }
end

local function channel_read(strand, opts)
	return channel_read_impl(strand, opts, false)
end

local function channel_read_soft(strand, opts)
	return channel_read_impl(strand, opts, true)
end

local function read_free(rd)
	if not rd then return end
	rd.strand = nil
	rd.qual = nil
	rd.bases = 0
	rd.dropped = 0
end

local function consensus_read(strand, opts)
	if opts == nil then return channel_read(strand, nil) end
	local ch = opts.ch
	local coverage = opts.coverage or 0
	local soft = opts.soft ~= nil and opts.soft ~= false and opts.soft ~= 0
	if coverage <= 1 then
		if soft then return channel_read_soft(strand, ch) end
		return channel_read(strand, ch)
	end
	if type(strand) ~= "string" then return nil, "expected string" end
	local n = #strand * 4
	if n == 0 then
		if soft then return channel_read_soft(strand, ch) end
		return channel_read(strand, ch)
	end
	local counts = {}
	for i = 0, n * 4 - 1 do counts[i] = 0 end
	local survivors = 0
	for r = 0, coverage - 1 do
		local c = {}
		if ch then for k, v in pairs(ch) do c[k] = v end end
		c.seed = ((ch and ch.seed or 0) + r) % TWO32
		local rd, err
		if soft then rd, err = channel_read_soft(strand, c) else rd, err = channel_read(strand, c) end
		if not rd then return nil, err end
		if rd.dropped ~= 1 then
			if rd.bases ~= n then
				read_free(rd)
				return nil, "consensus requires equal-length reads (no indels)"
			end
			for i = 0, n - 1 do
				local w = 1
				if soft and rd.qual then w = rd.qual[i + 1] + 1 end
				local idx = i * 4 + digit_at(rd.strand, i)
				counts[idx] = counts[idx] + w
			end
			survivors = survivors + 1
		end
		read_free(rd)
	end
	if survivors == 0 then
		return { strand = "", bases = 0, dropped = 1 }
	end
	local best = {}
	for i = 0, n - 1 do
		local b = 0
		for d = 1, 3 do
			if counts[i * 4 + d] > counts[i * 4 + b] then b = d end
		end
		best[i + 1] = b
	end
	return { strand = pack_digits(best, n), bases = n, dropped = 0 }
end

return {
	read = channel_read,
	read_soft = channel_read_soft,
	read_free = read_free,
	consensus_read = consensus_read,
	prng_new = prng_new,
	prng_next = prng_next,
	prng_next64 = prng_next64,
	prng_chance = prng_chance,
}
