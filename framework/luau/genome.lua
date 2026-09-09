--!native
-- genome.lua  (phase 1 -- the gene layer)
--
-- A living-DNA storage layer on top of dna.lua. Phase 1 implements the
-- GENE: a self-delimiting data segment framed by promoter/terminator
-- markers, so genes are found by scanning (no global length header),
-- like RNA polymerase reading a chromosome.
--
-- Gene layout (packed strands, 4 bases per byte):
--   [PROM: 12 bases "ACGTACGTACGT" = bytes 0x1B 0x1B 0x1B]
--   [header: 12 codons = 6 bytes: id, type, rawlen(16), packedlen(16)]
--   [payload: packedlen bytes]
--   [tag: 8 codons = 4 bytes, Chaskey-12 tag over the raw data]
--   [TERM: 12 bases "TGACTGACTGAC" = bytes 0xE4 0xE4 0xE4]
--
-- Codon code (mimics the degenerate genetic code):
--   * codon = 3 bases (b1, b2, wobble), 6 bits
--   * value = b1*4 + b2 (16 classes, 4 bits); the wobble base carries no
--     data -- like the third codon position in real DNA, any corruption
--     of it leaves the value unchanged (wobble is chosen != b2, which
--     also caps homopolymer runs at h = 3 inside codon sections)
--   * density: 4 bits per 3 bases (1.33 bits/nt); used for the header
--     and CRC (critical zones). Dense payloads are delegated to
--     dna.encode (~1.95 bits/nt) with the CRC guarding them
--   * type bit 0 selects the payload mode: 0 = dense, 1 = codon
--
-- Damage model: wobble corruption is silently tolerated (data intact,
-- CRC still valid); substitutions in b1/b2 or dense payload change the
-- data and are DETECTED by the CRC -- repairing them is phase 3
-- (diploid homologous recombination).
--
-- Requires bit32 (Luau, Lua 5.2+) and the dna module.

local dna = require("./dna")

local byte, char, concat, sub, find =
	string.byte, string.char, table.concat, string.sub, string.find
local floor = math.floor

local band, bxor, lshift, rshift, bor =
	bit32.band, bit32.bxor, bit32.lshift, bit32.rshift, bit32.bor

local PROM = "\27\27\27"    -- 12 bases: ACGT x3 (digit 0,1,2,3)
local TERM = "\228\228\228" -- 12 bases: TGAC x3 (digit 3,2,1,0)

-- Chaskey-12 (ISO/IEC 29192-6): the official MAC mode of Nicky Mouha's
-- reference implementation, with the tag truncated to 32 bits (the
-- design explicitly supports 32..128-bit tags). Unlike CRC32 there is
-- no linear structure, so ANY payload modification changes the tag with
-- probability ~1 - 2^-32 regardless of the corruption pattern. The
-- default key is fixed and public: this is an integrity tag, not a MAC
-- against active forgery.
local CHASKEY_KEY = { 0x1B873593, 0x9E3779B9, 0x85EBCA6B, 0xC2B2AE35 }
local C87 = { 0x00, 0x87 }

local function rotl(x, n)
	return bor(lshift(x, n), rshift(x, 32 - n))
end

-- CMAC-style subkey doubling over the 128-bit block (TimesTwo)
local function times_two(k)
	local c = rshift(k[4], 31)
	return {
		bxor(lshift(k[1], 1), C87[c + 1]),
		bor(lshift(k[2], 1), rshift(k[1], 31)),
		bor(lshift(k[3], 1), rshift(k[2], 31)),
		bor(lshift(k[4], 1), rshift(k[3], 31)),
	}
end

local function permute12(v0, v1, v2, v3)
	for _ = 1, 12 do
		v0 = band(v0 + v1, 0xFFFFFFFF)
		v1 = bxor(rotl(v1, 5), v0)
		v0 = rotl(v0, 16)
		v2 = band(v2 + v3, 0xFFFFFFFF)
		v3 = bxor(rotl(v3, 8), v2)
		v0 = band(v0 + v3, 0xFFFFFFFF)
		v3 = bxor(rotl(v3, 13), v0)
		v2 = band(v2 + v1, 0xFFFFFFFF)
		v1 = bxor(rotl(v1, 7), v2)
		v2 = rotl(v2, 16)
	end
	return v0, v1, v2, v3
end

local function load32(s, off)
	local b1, b2, b3, b4 = byte(s, off, off + 3)
	return bor(b1, lshift(b2, 8), lshift(b3, 16), lshift(b4, 24))
end

local function chaskey_tag(s, key)
	if type(s) ~= "string" then return nil, "expected string" end
	local k = key or CHASKEY_KEY
	local k1 = times_two(k)
	local k2 = times_two(k1)
	local len = #s
	local rem = len % 16
	local v0, v1, v2, v3 = k[1], k[2], k[3], k[4]
	local limit
	if rem == 0 then
		limit = len - 16  -- all full blocks except the last
	else
		limit = len - rem
	end
	for off = 1, limit, 16 do
		v0 = bxor(v0, load32(s, off))
		v1 = bxor(v1, load32(s, off + 4))
		v2 = bxor(v2, load32(s, off + 8))
		v3 = bxor(v3, load32(s, off + 12))
		v0, v1, v2, v3 = permute12(v0, v1, v2, v3)
	end
	local m0, m1, m2, m3
	local l0, l1, l2, l3
	if rem == 0 and len > 0 then
		local off = len - 15
		m0, m1, m2, m3 = load32(s, off), load32(s, off + 4),
			load32(s, off + 8), load32(s, off + 12)
		l0, l1, l2, l3 = k1[1], k1[2], k1[3], k1[4]
	else
		local p = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }
		for i = 1, rem do p[i] = byte(s, len - rem + i) end
		p[rem + 1] = 1  -- padding bit (10*)
		m0 = bor(p[1], lshift(p[2], 8), lshift(p[3], 16), lshift(p[4], 24))
		m1 = bor(p[5], lshift(p[6], 8), lshift(p[7], 16), lshift(p[8], 24))
		m2 = bor(p[9], lshift(p[10], 8), lshift(p[11], 16), lshift(p[12], 24))
		m3 = bor(p[13], lshift(p[14], 8), lshift(p[15], 16), lshift(p[16], 24))
		l0, l1, l2, l3 = k2[1], k2[2], k2[3], k2[4]
	end
	v0 = bxor(v0, m0, l0)
	v1 = bxor(v1, m1, l1)
	v2 = bxor(v2, m2, l2)
	v3 = bxor(v3, m3, l3)
	v0, v1, v2, v3 = permute12(v0, v1, v2, v3)
	return bxor(v0, l0), bxor(v1, l1), bxor(v2, l2), bxor(v3, l3)
end

local function chaskey32(s, key)
	local t0, err = chaskey_tag(s, key)
	if not t0 then return nil, err end
	return t0
end

-- codon packing: value = 4 bits -> codon (b1, b2, wobble)
local function pack_codons(values)
	local t = {}
	local acc, nib = 0, 0
	local function emit(d)
		acc = acc * 4 + d
		nib = nib + 1
		if nib == 4 then
			t[#t + 1] = char(acc)
			acc, nib = 0, 0
		end
	end
	for k = 1, #values do
		local c = values[k]
		local b1 = band(rshift(c, 2), 3)
		local b2 = band(c, 3)
		local w = bxor(b2, 1)  -- wobble: != b2, so every codon ends run 1
		emit(b1)
		emit(b2)
		emit(w)
	end
	if nib ~= 0 then return nil, "codon count must be a multiple of 4" end
	return concat(t)
end

local function read2(s, bitpos)
	-- 2-bit digit at bit offset `bitpos` of packed string s (offsets
	-- stay even for the 6-bit codon stride, but a straddle is handled
	-- for safety)
	local bytei = floor(bitpos / 8) + 1
	local off = bitpos % 8
	local b = byte(s, bytei)
	if off <= 6 then
		return band(rshift(b, 6 - off), 3)
	end
	local nb = byte(s, bytei + 1) or 0
	return bor(lshift(band(b, 1), 1), rshift(nb, 7))
end

local function unpack_codons(packed, nvalues)
	local vals = {}
	local bitpos = 0
	for k = 1, nvalues do
		vals[k] = read2(packed, bitpos) * 4 + read2(packed, bitpos + 2)
		bitpos = bitpos + 6
	end
	return vals
end

local function values_from_bytes(s)
	local v = {}
	for i = 1, #s do
		local b = byte(s, i)
		v[#v + 1] = band(rshift(b, 4), 15)
		v[#v + 1] = band(b, 15)
	end
	return v
end

local function bytes_from_values(v, nbytes)
	local t = {}
	for i = 1, nbytes do
		t[#t + 1] = char(v[2 * i - 1] * 16 + v[2 * i])
	end
	return concat(t)
end

local function gene_encode(data, opts)
	opts = opts or {}
	if type(data) ~= "string" then return nil, "expected string" end
	local id = opts.id or 0
	local usertype = opts.type or 0
	local mode = opts.mode or "dense"
	if type(id) ~= "number" or id ~= floor(id) or id < 0 or id > 255 then
		return nil, "invalid id (byte 0..255)"
	end
	if type(usertype) ~= "number" or usertype ~= floor(usertype)
		or usertype < 0 or usertype > 255 then
		return nil, "invalid type (byte 0..255)"
	end
	if mode ~= "dense" and mode ~= "codon" then
		return nil, "invalid mode (dense|codon)"
	end
	local rawlen = #data
	if rawlen > 65535 then return nil, "gene too large (raw > 65535)" end
	local payload, packedlen
	if mode == "dense" then
		local h = opts.h or 3
		if type(h) ~= "number" or h ~= floor(h) or h < 3 or h > 12 then
			return nil, "invalid h for gene (integer 3..12)"
		end
		payload = dna.encode(data, { h = h })
		if not payload then return nil, "dense payload encode failed" end
		packedlen = #payload
	else
		local v = values_from_bytes(data)
		while #v % 4 ~= 0 do v[#v + 1] = 0 end
		payload = pack_codons(v)
		packedlen = #payload
	end
	if packedlen > 65535 then return nil, "gene too large (packed > 65535)" end
	local typ = bor(band(usertype, 0xFE), (mode == "codon") and 1 or 0)
	local hdr = char(id, typ,
		band(rshift(rawlen, 8), 255), band(rawlen, 255),
		band(rshift(packedlen, 8), 255), band(packedlen, 255))
	local hdrp = pack_codons(values_from_bytes(hdr))
	local c = chaskey32(data)
	if c == nil then return nil, "tag computation failed" end
	local crcp = pack_codons(values_from_bytes(char(
		band(rshift(c, 24), 255), band(rshift(c, 16), 255),
		band(rshift(c, 8), 255), band(c, 255))))
	return PROM .. hdrp .. payload .. crcp .. TERM
end

local function parse_gene_at(strand, p)
	if p + 20 > #strand then return nil end
	local hv = unpack_codons(sub(strand, p + 3, p + 11), 12)
	-- header values are 4-bit nibbles: two per byte
	local id = hv[1] * 16 + hv[2]
	local typ = hv[3] * 16 + hv[4]
	local rawlen = hv[5] * 4096 + hv[6] * 256 + hv[7] * 16 + hv[8]
	local packedlen = hv[9] * 4096 + hv[10] * 256 + hv[11] * 16 + hv[12]
	if p + 20 + packedlen > #strand then return nil end
	local termi = p + 18 + packedlen
	if sub(strand, termi, termi + 2) ~= TERM then return nil end
	local modebit = band(typ, 1)
	local payload = sub(strand, p + 12, p + 11 + packedlen)
	local cv = unpack_codons(sub(strand, p + 12 + packedlen, p + 17 + packedlen), 8)
	local expected = (cv[1] * 16 + cv[2]) * 0x1000000
		+ (cv[3] * 16 + cv[4]) * 0x10000
		+ (cv[5] * 16 + cv[6]) * 0x100
		+ (cv[7] * 16 + cv[8])
	local data, crc_ok = nil, false
	if modebit == 0 then
		local raw = dna.decode(payload)
		if raw then
			data = sub(raw, 1, rawlen)
			crc_ok = (#data == rawlen) and chaskey32(data) == expected
		end
	else
		local nv = packedlen * 4 / 3
		if nv % 4 == 0 and nv % 1 == 0 then
			local raw = bytes_from_values(unpack_codons(payload, nv), nv / 2)
			data = sub(raw, 1, rawlen)
			crc_ok = (#data == rawlen) and chaskey32(data) == expected
		end
	end
	return {
		id = id,
		type = typ,
		mode = (modebit == 1) and "codon" or "dense",
		rawlen = rawlen,
		packedlen = packedlen,
		offset = p,
		size = 21 + packedlen,
		data = data,
		crc_ok = crc_ok,
		tag32 = expected,
	}
end

local function scan(strand)
	if type(strand) ~= "string" then return nil, "expected string" end
	local genes = {}
	local pos = 1
	while pos <= #strand do
		local p = find(strand, PROM, pos, true)
		if not p then break end
		local g = parse_gene_at(strand, p)
		if g then
			genes[#genes + 1] = g
			pos = p + g.size
		else
			pos = p + 1
		end
	end
	return genes
end

local function gene_read(strand, id)
	local genes, err = scan(strand)
	if not genes then return nil, err end
	for i = 1, #genes do
		local g = genes[i]
		if g.crc_ok and (id == nil or g.id == id) then
			return g.data, g
		end
	end
	if id ~= nil then return nil, "no valid gene with id " .. id end
	return nil, "no valid gene"
end

-- rewrite a single base (0-based index); used by damage simulators
local function set_base(strand, idx, d)
	local bytei = floor(idx / 4) + 1
	local sh = 2 * (3 - idx % 4)
	local old = byte(strand, bytei)
	local nb = bor(band(old, bxor(255, lshift(3, sh))), lshift(band(d, 3), sh))
	return sub(strand, 1, bytei - 1) .. char(nb) .. sub(strand, bytei + 1)
end

return {
	chaskey32 = chaskey32,
	chaskey_tag = chaskey_tag,
	pack_codons = pack_codons,
	unpack_codons = unpack_codons,
	values_from_bytes = values_from_bytes,
	bytes_from_values = bytes_from_values,
	gene_encode = gene_encode,
	scan = scan,
	gene_read = gene_read,
	set_base = set_base,
	PROM = PROM,
	TERM = TERM,
}

