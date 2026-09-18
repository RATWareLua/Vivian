--!native
local byte, char, concat = string.byte, string.char, table.concat
local floor = math.floor

local band, bxor, lshift, rshift = bit32.band, bit32.bxor, bit32.lshift, bit32.rshift

local GF_POLY = 0x11D
local RS_MAX_N = 255

local gf_exp, gf_log = {}, {}
do
	local x = 1
	for i = 0, 254 do
		gf_exp[i] = x
		gf_log[x] = i
		x = lshift(x, 1)
		if band(x, 0x100) ~= 0 then x = bxor(x, GF_POLY) end
	end
	for i = 255, 511 do gf_exp[i] = gf_exp[i - 255] end
	gf_log[0] = 0
end

local function gf_mul(a, b)
	if a == 0 or b == 0 then return 0 end
	local t = gf_log[a] + gf_log[b]
	if t >= 255 then t = t - 255 end
	return gf_exp[t]
end

local function gf_div(a, b)
	local t = gf_log[a] - gf_log[b]
	if t < 0 then t = t + 255 end
	return gf_exp[t]
end

local function gf_inv(a)
	return gf_exp[255 - gf_log[a]]
end

local function gf_pow(a, e)
	if a == 0 then return 0 end
	return gf_exp[(gf_log[a] * e) % 255]
end

local function poly_eval(v, n, x)
	local y = 0
	for i = 0, n - 1 do y = bxor(gf_mul(y, x), v[i]) end
	return y
end

local function poly_mul(a, na, b, nb, out)
	local nc = na + nb - 1
	for i = 0, nc - 1 do out[i] = 0 end
	for i = 0, na - 1 do
		for j = 0, nb - 1 do
			out[i + j] = bxor(out[i + j], gf_mul(a[i], b[j]))
		end
	end
end

local function lagrange_interp(xs, ys, m, out)
	local basis, tmp = {}, {}
	for i = 0, m - 1 do out[i] = 0 end
	for j = 0, m - 1 do
		local nb = 1
		basis[0] = 1
		local den = 1
		for l = 0, m - 1 do
			if l ~= j then
				local lin = { [0] = xs[l], [1] = 1 }
				poly_mul(basis, nb, lin, 2, tmp)
				nb = nb + 1
				for i = 0, nb - 1 do basis[i] = tmp[i] end
				den = gf_mul(den, bxor(xs[j], xs[l]))
			end
		end
		local scale = gf_mul(ys[j], gf_inv(den))
		for t = 0, m - 1 do out[t] = bxor(out[t], gf_mul(basis[t], scale)) end
	end
end

local function encode(data, k, m)
	if data ~= nil and type(data) ~= "string" then return nil, "expected string" end
	if data == nil and k ~= 0 then return nil, "expected buffers" end
	if m < 1 or k + m > RS_MAX_N then return nil, "invalid rs geometry" end
	local d = {}
	if k > 0 then
		for i = 0, k - 1 do d[i] = byte(data, i + 1) end
	end
	local xs, ys, parity = {}, {}, {}
	for j = 0, m - 1 do
		xs[j] = gf_exp[j]
		ys[j] = gf_mul(poly_eval(d, k, xs[j]), gf_pow(xs[j], m))
	end
	lagrange_interp(xs, ys, m, parity)
	local t = {}
	if k > 0 then t[#t + 1] = data end
	for r = 0, m - 1 do t[#t + 1] = char(parity[m - 1 - r]) end
	return concat(t)
end

local function syndromes(cw, n, m, S)
	for j = 0, m - 1 do S[j] = poly_eval(cw, n, gf_exp[j]) end
end

local function berlekamp_massey(S, m, C)
	local B, T = {}, {}
	for i = 0, m - 1 do C[i] = 0; B[i] = 0 end
	C[0] = 1
	B[0] = 1
	local L, shift, b = 0, 1, 1
	for nn = 0, m - 1 do
		local d = S[nn]
		local lim = (L < nn) and L or nn
		for i = 1, lim do d = bxor(d, gf_mul(C[i], S[nn - i])) end
		if d == 0 then
			shift = shift + 1
		else
			local coef = gf_div(d, b)
			if 2 * L <= nn then
				for i = 0, m - 1 do T[i] = C[i] end
				for i = 0, m - 1 - shift do
					C[i + shift] = bxor(C[i + shift], gf_mul(coef, B[i]))
				end
				L = nn + 1 - L
				for i = 0, m - 1 do B[i] = T[i] end
				b = d
				shift = 1
			else
				for i = 0, m - 1 - shift do
					C[i + shift] = bxor(C[i + shift], gf_mul(coef, B[i]))
				end
				shift = shift + 1
			end
		end
	end
	return L
end

local function solve_values(S, X, L, e)
	local A = {}
	for j = 0, L - 1 do A[j] = {} end
	for l = 0, L - 1 do
		local pw = 1
		for j = 0, L - 1 do
			A[j][l] = pw
			pw = gf_mul(pw, X[l])
		end
	end
	for j = 0, L - 1 do A[j][L] = S[j] end
	for col = 0, L - 1 do
		local piv = L
		for r = col, L - 1 do
			if A[r][col] ~= 0 then piv = r; break end
		end
		if piv == L then return false end
		if piv ~= col then
			for c = col, L do
				local t = A[col][c]; A[col][c] = A[piv][c]; A[piv][c] = t
			end
		end
		local inv = gf_inv(A[col][col])
		for c = col, L do A[col][c] = gf_mul(A[col][c], inv) end
		for r = 0, L - 1 do
			if r ~= col and A[r][col] ~= 0 then
				local f = A[r][col]
				for c = col, L do
					A[r][c] = bxor(A[r][c], gf_mul(f, A[col][c]))
				end
			end
		end
	end
	for i = 0, L - 1 do e[i] = A[i][L] end
	return true
end

local function decode(cw, n, k)
	if type(cw) ~= "string" then return nil, "expected string" end
	if n > RS_MAX_N or k < 1 or k >= n then return nil, "invalid rs geometry" end
	local m = n - k
	local bytes = {}
	for i = 0, n - 1 do bytes[i] = byte(cw, i + 1) end
	local S = {}
	syndromes(bytes, n, m, S)
	local nonzero = false
	for j = 0, m - 1 do
		if S[j] ~= 0 then nonzero = true; break end
	end
	if not nonzero then return cw, 0 end
	local C = {}
	local L = berlekamp_massey(S, m, C)
	if L == 0 or L > floor(m / 2) then return nil, "too many errors" end
	local pos, X, e = {}, {}, {}
	local found = 0
	for i = 0, n - 1 do
		local x = gf_exp[n - 1 - i]
		local xinv = gf_inv(x)
		local val, pw = 0, 1
		for t = 0, L do
			val = bxor(val, gf_mul(C[t], pw))
			pw = gf_mul(pw, xinv)
		end
		if val == 0 then
			if found >= L then return nil, "too many roots" end
			pos[found] = i
			X[found] = x
			found = found + 1
		end
	end
	if found ~= L then return nil, "error locator has no roots" end
	if not solve_values(S, X, L, e) then return nil, "singular system" end
	for l = 0, L - 1 do bytes[pos[l]] = bxor(bytes[pos[l]], e[l]) end
	local S2 = {}
	syndromes(bytes, n, m, S2)
	for j = 0, m - 1 do
		if S2[j] ~= 0 then return nil, "miscorrection" end
	end
	local out = {}
	for i = 0, n - 1 do out[#out + 1] = char(bytes[i]) end
	return concat(out), L
end

return {
	encode = encode,
	decode = decode,
	gf_mul = gf_mul,
	gf_div = gf_div,
	gf_inv = gf_inv,
	gf_pow = gf_pow,
	poly_eval = poly_eval,
	poly_mul = poly_mul,
	lagrange_interp = lagrange_interp,
}
