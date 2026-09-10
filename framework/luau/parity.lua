--!native
-- parity.lua  (research layer -- systematic Reed-Solomon over GF(256))
--
-- Byte-identical port of src/parity.c: Lagrange-form systematic RS with
-- the classic 0x11D polynomial and generator 2. Shards are Lua strings,
-- so encoding and decoding are pure and side-effect free.
--
-- encode(data, m, len) -> n+m shards (data shards copied 1:1)
-- decode(shards, present, n, m, len) -> n data shards
--
-- Requires bit32 (Luau, Lua 5.2+).

local byte, char, concat, sub, unpack =
	string.byte, string.char, table.concat, string.sub, table.unpack or unpack
local band, bxor, lshift = bit32.band, bit32.bxor, bit32.lshift

local GF_POLY = 0x11D

local gf_exp, gf_log

local function gf_init()
	if gf_exp then return end
	gf_exp, gf_log = {}, {}
	local x = 1
	for i = 0, 254 do
		gf_exp[i] = x
		gf_log[x] = i
		x = lshift(x, 1)
		if band(x, 0x100) ~= 0 then x = bxor(x, GF_POLY) end
	end
	for i = 255, 511 do gf_exp[i] = gf_exp[i - 255] end
end

local function gf_mul(a, b)
	if a == 0 or b == 0 then return 0 end
	local t = gf_log[a] + gf_log[b]
	if t >= 255 then t = t - 255 end
	return gf_exp[t]
end

local function gf_inv(a)
	return gf_exp[255 - gf_log[a]]
end

-- Lagrange basis coefficients at x over the points xs[1..n]
local function lagrange_coefs(xs, n, x)
	local coef = {}
	for i = 1, n do
		local num, den = 1, 1
		for k = 1, n do
			if k ~= i then
				num = gf_mul(num, bxor(x, xs[k]))
				den = gf_mul(den, bxor(xs[i], xs[k]))
			end
		end
		coef[i] = gf_mul(num, gf_inv(den))
	end
	return coef
end

local function bytes_to_string(t)
	local res = {}
	for i = 1, #t, 2048 do
		local j = (i + 2047 <= #t) and (i + 2047) or #t
		res[#res + 1] = char(unpack(t, i, j))
	end
	return concat(res)
end

-- data[1..n] are the data shards (each at least len bytes); returns
-- n+m shards: the data copies followed by m parity shards
local function encode(data, m, len)
	gf_init()
	local n = #data
	if n < 1 or len < 1 or n + m > 255 then
		return nil, "invalid shard geometry"
	end
	local xs = {}
	for i = 1, n do xs[i] = i end
	local shards = {}
	for i = 1, n do
		local s = data[i]
		if type(s) ~= "string" or #s < len then return nil, "invalid shard" end
		shards[i] = sub(s, 1, len)
	end
	for j = 1, m do
		local coef = lagrange_coefs(xs, n, n + j)
		local t = {}
		for b = 1, len do
			local acc = 0
			for i = 1, n do
				acc = bxor(acc, gf_mul(coef[i], byte(data[i], b)))
			end
			t[b] = acc
		end
		shards[n + j] = bytes_to_string(t)
	end
	return shards
end

-- shards[1..n+m] as strings (nil = erased), present[1..n+m] flags;
-- returns the n data shards
local function decode(shards, present, n, m, len)
	gf_init()
	if n < 1 or len < 1 or n + m > 255 then
		return nil, "invalid shard geometry"
	end
	local px = {}
	local np = 0
	for k = 1, n + m do
		if present[k] then
			if not shards[k] then return nil, "present shard is nil" end
			np = np + 1
			px[np] = k
			if np >= n then break end
		end
	end
	if np < n then return nil, "too few shards" end
	local out = {}
	for i = 1, n do
		if present[i] then
			out[i] = sub(shards[i], 1, len)
		else
			local coef = lagrange_coefs(px, n, i)
			local t = {}
			for b = 1, len do
				local acc = 0
				for k = 1, n do
					acc = bxor(acc, gf_mul(coef[k], byte(shards[px[k]], b)))
				end
				t[b] = acc
			end
			out[i] = bytes_to_string(t)
		end
	end
	return out
end

return {
	encode = encode,
	decode = decode,
}
