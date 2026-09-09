--!native
-- dna.lua  (v2 -- packed quaternary codec)
--
-- Quaternary (DNA-style) in-memory storage codec. Data is stored as a
-- packed strand: one nucleotide per 2 bits, 4 bases per byte, over the
-- alphabet A/C/G/T (00 = A, 01 = C, 10 = G, 11 = T), so a strand costs
-- ~2% more memory than the raw payload. ASCII views ("ACGT" strings)
-- are available via to_ascii/from_ascii.

local byte, char, concat, sub =
	string.byte, string.char, table.concat, string.sub
local floor = math.floor
local unpack = table.unpack or unpack

local band, bxor, lshift, rshift, bor =
	bit32.band, bit32.bxor, bit32.lshift, bit32.rshift, bit32.bor

local CHARS = { "A", "C", "G", "T" }  -- digits 0..3
local BDIGITS = { [65] = 0, [67] = 1, [71] = 2, [84] = 3 }  -- A C G T

-- ASCII[b]: 4-char ACGT view of packed byte b
local ASCII = {}
for b = 0, 255 do
	ASCII[b] = CHARS[band(rshift(b, 6), 3) + 1]
		.. CHARS[band(rshift(b, 4), 3) + 1]
		.. CHARS[band(rshift(b, 2), 3) + 1]
		.. CHARS[band(b, 3) + 1]
end

local function xs32(x)
	x = bxor(x, lshift(x, 13))
	x = bxor(x, rshift(x, 17))
	x = bxor(x, lshift(x, 5))
	return x
end

local function seed_for(len)
	local s = bxor(bxor(0x9E3779B9, len), 0x5BD1E995)
	-- xorshift32 is stuck at zero; length 0xC5E6902C would seed exactly 0
	return xs32((s == 0) and 1 or s)
end

local HDR_SEED = xs32(0x1B873593)

-- keystream bytes: each byte is the next 8 keystream bits, MSB first
local function make_keystream(seed)
	local state = seed
	local word, kleft = 0, 0
	return function()
		if kleft == 0 then
			state = xs32(state)
			word = state
			kleft = 32
		end
		kleft = kleft - 8
		return band(rshift(word, kleft), 255)
	end
end

-- PAIRS[boff * 256 + b]: the 2-bit digit starting at bit offset `boff`
-- of byte b (offsets 0..6; a 2-bit read never starts at offset 7)
local PAIRS = {}
for b = 0, 255 do
	for boff = 0, 6 do
		PAIRS[boff * 256 + b] = band(rshift(b, 6 - boff), 3)
	end
end

-- GCT[b]: number of C/G digits (1, 2) packed in byte b
local GCT = {}
for b = 0, 255 do
	local g, x = 0, b
	for _ = 1, 4 do
		local d = band(x, 3)
		if d == 1 or d == 2 then g = g + 1 end
		x = rshift(x, 2)
	end
	GCT[b] = g
end

-- BINFO[b] = { sd, sr, mr, ed, er }: start digit, start run, max run,
-- end digit and end run of the 4 digits packed in byte b
local BINFO = {}
for b = 0, 255 do
	local d0 = band(rshift(b, 6), 3)
	local d1 = band(rshift(b, 4), 3)
	local d2 = band(rshift(b, 2), 3)
	local d3 = band(b, 3)
	local mr, run, prev = 1, 1, d0
	if d1 == prev then run = run + 1 else run, prev = 1, d1 end
	if run > mr then mr = run end
	if d2 == prev then run = run + 1 else run, prev = 1, d2 end
	if run > mr then mr = run end
	if d3 == prev then run = run + 1 else run, prev = 1, d3 end
	if run > mr then mr = run end
	local sr = 1
	if d1 == d0 then
		sr = 2
		if d2 == d0 then
			sr = 3
			if d3 == d0 then sr = 4 end
		end
	end
	local er = 1
	if d3 == d2 then
		er = 2
		if d2 == d1 then
			er = 3
			if d1 == d0 then er = 4 end
		end
	end
	BINFO[b] = { sd = d0, sr = sr, mr = mr, ed = d3, er = er }
end

-- Byte fast path: when all 4 digits of a byte decode in 2-bit mode (no
-- run reaches h at a step start), the state machine collapses to one
-- table entry. Value layout: 1 + nrun*128 + nlast*16 + ngc*2, where
-- nrun/nlast is the state after the byte and ngc the digits that are C/G
local FASTCACHE = {}

local function fast_table(h)
	local t = FASTCACHE[h]
	if t then return t end
	t = {}
	for run = 1, h - 1 do
		for last = 0, 3 do
			local base = ((run - 1) * 4 + last) * 256
			for b = 0, 255 do
				local r, p, g, good = run, last, 0, true
				for sh = 6, 0, -2 do
					local d = band(rshift(b, sh), 3)
					if r >= h then
						good = false
						break
					end
					if d == 1 or d == 2 then g = g + 1 end
					if d == p then r = r + 1 else r, p = 1, d end
				end
				if good then
					t[base + b] = 1 + r * 128 + p * 16 + g * 2
				end
			end
		end
	end
	FASTCACHE[h] = t
	return t
end

local function encode(data, opts)
	if type(data) ~= "string" then return nil, "expected string" end
	if opts ~= nil and type(opts) ~= "table" then return nil, "expected table for opts" end
	local h = opts and opts.h or 3
	local eps = opts and opts.gc_eps or 0.05
	if type(h) ~= "number" or h ~= floor(h) or h < 1 or h > 12 then
		return nil, "invalid h (integer 1..12)"
	end
	if type(eps) ~= "number" or eps < 0.005 or eps >= 0.5 then
		return nil, "invalid gc_eps ([0.005, 0.5))"
	end
	local len = #data
	if len > 0xFFFFFFFF then return nil, "payload exceeds 4GB limit" end
	local nbits = 32 + 8 * len
	local hp = h - 1
	local d0 = hp % 4
	local d1 = (floor(hp / 4) + d0 + 1) % 4
	local kh = make_keystream(HDR_SEED)
	local kp = make_keystream(seed_for(len))

	local WHB = {
		bxor(band(rshift(len, 24), 255), kh()),
		bxor(band(rshift(len, 16), 255), kh()),
		bxor(band(rshift(len, 8), 255), kh()),
		bxor(band(len, 255), kh()),
	}
	local pj, pbc = 0, 0
	local function wfetch(w)
		if w <= 4 then return WHB[w] end
		local j = w - 4
		if j > len then return 0 end
		if j ~= pj then
			pj = j
			pbc = bxor(byte(data, j), kp())
		end
		return pbc
	end
	local out = {}
	local FT = fast_table(h)
	local acc, nib = d0 * 4 + d1, 2
	local gc = ((d0 == 1 or d0 == 2) and 1 or 0) + ((d1 == 1 or d1 == 2) and 1 or 0)
	local wi, wb, boff = 1, WHB[1], 0
	local kidx = wb
	local last, run = d1, h
	local bi = 0
	while bi < nbits do
		local remaining = (bi < 32) and (32 - bi) or (nbits - bi)
		local d
		local fast = false
		if run < h and boff == 0 and nib == 0 and (bi >= 32 or bi + 8 <= 32) then
			local e = FT[((run - 1) * 4 + last) * 256 + wb]
			if e then
				out[#out + 1] = char(wb)
				local q = e - 1
				run = band(rshift(q, 7), 15)
				last = band(rshift(q, 4), 3)
				gc = gc + band(rshift(q, 1), 7)
				bi = bi + 8
				wi = wi + 1
				wb = wfetch(wi)
				kidx = wb
				fast = true
			end
		end
		if not fast then
			if run >= h or remaining == 1 then
				local bit = band(rshift(wb, 7 - boff), 1)
				boff = boff + 1
				if boff == 8 then
					wi = wi + 1
					wb = wfetch(wi)
					boff = 0
				end
				if last == 0 or last == 3 then
					d = (bit == 0) and 1 or 2
				else
					d = (bit == 0) and 0 or 3
				end
				run, last = 1, d
				bi = bi + 1
				kidx = boff * 256 + wb
			else
				if boff <= 6 then
					d = PAIRS[kidx]
					kidx = kidx + 512
					if kidx >= 2048 then
						wi = wi + 1
						wb = wfetch(wi)
						kidx = wb
					end
					boff = boff + 2
					if boff == 8 then boff = 0 end
				else
					local nb = wfetch(wi + 1)
					d = bor(lshift(band(wb, 1), 1), rshift(nb, 7))
					wi = wi + 1
					wb = nb
					boff = 1
					kidx = 256 + nb
				end
				bi = bi + 2
				if d == last then run = run + 1 else run, last = 1, d end
			end
			acc = acc * 4 + d
			nib = nib + 1
			if nib == 4 then
				out[#out + 1] = char(acc)
				acc, nib = 0, 0
			end
			if d == 1 or d == 2 then gc = gc + 1 end
		end
	end

	local n = #out * 4 + nib
	local lo_b, hi_b = 0.5 - eps, 0.5 + eps
	local AT, GC = { 0, 3 }, { 1, 2 }
	local function pad(count, pair)
		for _ = 1, count do
			local d = (last == pair[1]) and pair[2] or pair[1]
			acc = acc * 4 + d
			nib = nib + 1
			if nib == 4 then
				out[#out + 1] = char(acc)
				acc, nib = 0, 0
			end
			n, last = n + 1, d
			if d == 1 or d == 2 then gc = gc + 1 end
		end
	end
	while true do
		local nf = (nib > 0) and (4 - nib) or 0
		local done = false
		local cands = (nf > 0) and { AT, GC } or { AT }
		for ci2 = 1, #cands do
			local pair = cands[ci2]
			local r = (gc + (pair == GC and nf or 0)) / (n + nf)
			if r >= lo_b and r <= hi_b then
				if nf > 0 then pad(nf, pair) end
				done = true
				break
			end
		end
		if done then break end
		local below = gc / n < 0.5
		local k = (nf > 0) and nf or 4
		local r2 = (gc + (below and k or 0)) / (n + k)
		if below and r2 > hi_b then
			pad(1, GC)
		elseif not below and r2 < lo_b then
			pad(1, AT)
		else
			pad(k, below and GC or AT)
		end
	end
	return concat(out)
end

local function decode(packed)
	if type(packed) ~= "string" then return nil, "expected string" end
	local nraw = #packed
	if nraw < 1 then return nil, "invalid strand" end
	local b0 = byte(packed, 1)
	local d0 = band(rshift(b0, 6), 3)
	local d1 = band(rshift(b0, 4), 3)
	local fc = (d1 - d0 - 1) % 4
	if fc > 2 then return nil, "invalid strand prefix" end
	local h = d0 + 4 * fc + 1
	local kh = make_keystream(HDR_SEED)
	local kp
	local bytes = {}
	local need_bytes
	local acc8, nacc = 0, 0
	local nbitsout = 0
	local function flush(w)
		bytes[#bytes + 1] = bxor(w, (#bytes < 4) and kh() or kp())
		acc8, nacc = 0, 0
		if not need_bytes and #bytes == 4 then
			need_bytes = 4 + bytes[1] * 0x1000000
				+ bytes[2] * 0x10000 + bytes[3] * 0x100 + bytes[4]
			kp = make_keystream(seed_for(need_bytes - 4))
		end
	end
	local FT = fast_table(h)
	local cb, ci = b0, 1
	local kidx = 1024 + b0
	local total_bases = nraw * 4
	local k = 2
	local last, run = d1, h
	while true do
		if need_bytes and #bytes >= need_bytes then break end
		if k >= total_bases then return nil, "truncated strand" end
		local before = nbitsout
		local need_bits = (need_bytes and need_bytes * 8) or 32
		local remaining = (before < 32) and (32 - before) or (need_bits - before)
		local fast = false
		if run < h and nacc == 0 and kidx == cb and remaining >= 8
			and (before >= 32 or before + 8 <= 32) then
			local e = FT[((run - 1) * 4 + last) * 256 + cb]
			if e then
				flush(cb)
				local q = e - 1
				run = band(rshift(q, 7), 15)
				last = band(rshift(q, 4), 3)
				nbitsout = before + 8
				k = k + 4
				ci = ci + 1
				cb = (ci <= nraw) and byte(packed, ci) or 0
				kidx = cb
				fast = true
			end
		end
		if not fast then
			local d = PAIRS[kidx]
			kidx = kidx + 512
			if kidx >= 2048 then
				ci = ci + 1
				cb = (ci <= nraw) and byte(packed, ci) or 0
				kidx = cb
			end
			k = k + 1
			before = nbitsout
			need_bits = (need_bytes and need_bytes * 8) or 32
			remaining = (before < 32) and (32 - before) or (need_bits - before)
			if run >= h or remaining == 1 then
				if d == last then return nil, "homopolymer violation" end
				local bt
				if last == 0 or last == 3 then
					if d ~= 1 and d ~= 2 then return nil, "invalid transition" end
					bt = d - 1
				else
					if d ~= 0 and d ~= 3 then return nil, "invalid transition" end
					bt = (d == 3) and 1 or 0
				end
				acc8 = acc8 * 2 + bt
				nacc = nacc + 1
				nbitsout = before + 1
				run, last = 1, d
				if nacc == 8 then flush(acc8) end
			else
				if nacc == 7 then
					acc8 = acc8 * 2 + band(rshift(d, 1), 1)
					flush(acc8)
					acc8 = band(d, 1)
					nacc = 1
				else
					acc8 = acc8 * 4 + d
					nacc = nacc + 2
					if nacc == 8 then flush(acc8) end
				end
				nbitsout = before + 2
				if d == last then run = run + 1 else run, last = 1, d end
			end
		end
	end

	-- Fast chunked string assembly avoids millions of 1-char allocations
	local res = {}
	for i = 5, #bytes, 2048 do
		local j = (i + 2047 <= #bytes) and (i + 2047) or #bytes
		res[#res + 1] = char(unpack(bytes, i, j))
	end
	return concat(res)
end

local function to_ascii(packed, nbases)
	if type(packed) ~= "string" then return nil, "expected string" end
	if nbases ~= nil and (type(nbases) ~= "number" or nbases ~= floor(nbases) or nbases < 0) then
		return nil, "invalid nbases (non-negative integer)"
	end
	nbases = nbases or #packed * 4
	local full = floor(nbases / 4)
	if full > #packed then full = #packed end
	local t = {}
	for i = 1, full do
		t[#t + 1] = ASCII[byte(packed, i)]
	end
	local r = nbases - full * 4
	if r > 0 and r < 4 and full < #packed then
		t[#t + 1] = sub(ASCII[byte(packed, full + 1)], 1, r)
	end
	return concat(t)
end

local function from_ascii(s)
	if type(s) ~= "string" then return nil, "expected string" end
	s = (s:gsub("%s+", "")):upper()
	if s:find("[^ACGT]") then return nil, "invalid base" end
	local t, acc, k = {}, 0, 0
	for i = 1, #s do
		acc = acc * 4 + BDIGITS[byte(s, i)]
		k = k + 1
		if k == 4 then
			t[#t + 1] = char(acc)
			acc, k = 0, 0
		end
	end
	if k > 0 then t[#t + 1] = char(lshift(acc, (4 - k) * 2)) end
	return concat(t)
end

local REV = {}
for b = 0, 255 do
	REV[b] = bor(
		lshift(band(b, 3), 6),
		lshift(band(rshift(b, 2), 3), 4),
		lshift(band(rshift(b, 4), 3), 2),
		band(rshift(b, 6), 3))
end

local function complement(packed)
	if type(packed) ~= "string" then return nil, "expected string" end
	local t = {}
	for i = 1, #packed do t[i] = char(bxor(byte(packed, i), 0xFF)) end
	return concat(t)
end

local function reverse_complement(packed)
	if type(packed) ~= "string" then return nil, "expected string" end
	local n = #packed
	local t = {}
	for i = 1, n do t[i] = char(bxor(REV[byte(packed, n + 1 - i)], 0xFF)) end
	return concat(t)
end

local function gc_content(packed)
	if type(packed) ~= "string" then return nil, "expected string" end
	local n = #packed
	if n == 0 then return 0 end
	local gc = 0
	for i = 1, n do gc = gc + GCT[byte(packed, i)] end
	return gc / (n * 4)
end

local function max_homopolymer(packed)
	if type(packed) ~= "string" then return nil, "expected string" end
	local best, prevd, prevr = 0, -1, 0
	for i = 1, #packed do
		local inf = BINFO[byte(packed, i)]
		if inf.mr > best then best = inf.mr end
		if inf.sd == prevd and prevd >= 0 then
			local r = prevr + inf.sr
			if r > best then best = r end
			prevr = (inf.sr == 4) and (prevr + 4) or inf.er
		else
			prevr = inf.er
		end
		prevd = inf.ed
	end
	return best
end

local function validate(packed, opts)
	if type(packed) ~= "string" then return nil, "expected string" end
	if #packed < 1 then return nil, "invalid strand" end
	if opts ~= nil and type(opts) ~= "table" then return nil, "expected table for opts" end
	local eps = opts and opts.gc_eps or 0.05
	if type(eps) ~= "number" or eps < 0.005 or eps >= 0.5 then
		return nil, "invalid gc_eps ([0.005, 0.5))"
	end
	local h
	if opts and opts.h ~= nil then
		h = opts.h
		if type(h) ~= "number" or h ~= floor(h) or h < 1 or h > 12 then
			return nil, "invalid h (integer 1..12)"
		end
	else
		local b0 = byte(packed, 1)
		local d0 = band(rshift(b0, 6), 3)
		local d1 = band(rshift(b0, 4), 3)
		local fc = (d1 - d0 - 1) % 4
		if fc > 2 then return nil, "invalid strand prefix" end
		h = d0 + 4 * fc + 1
	end
	local m = max_homopolymer(packed)
	if m > h then
		return nil, "homopolymer run " .. m .. " exceeds h = " .. h
	end
	local g = gc_content(packed)
	if g < 0.5 - eps or g > 0.5 + eps then
		return nil, string.format(
			"gc content %.3f outside [%.3f, %.3f]", g, 0.5 - eps, 0.5 + eps)
	end
	return true
end

return {
	encode = encode,
	decode = decode,
	to_ascii = to_ascii,
	from_ascii = from_ascii,
	complement = complement,
	reverse_complement = reverse_complement,
	gc_content = gc_content,
	max_homopolymer = max_homopolymer,
	validate = validate,
}