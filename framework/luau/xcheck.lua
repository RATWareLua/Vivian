--!native
-- xcheck.lua -- deterministic Lua side of the C<->Lua byte-parity check.
--
-- Prints one labeled hex line per scenario; test/xcheck.c prints the
-- same lines and CI diffs the two outputs. Keep the scenario list,
-- payload seeds and all option values in sync with the C file.

local chromosome = require("./chromosome")
local organism = require("./organism")

local byte, char, concat, sub =
	string.byte, string.char, table.concat, string.sub
local band, bxor, rshift, unpack =
	bit32.band, bit32.bxor, bit32.rshift, table.unpack or unpack

local function bytes_to_string(t)
	local res = {}
	for i = 1, #t, 2048 do
		local j = (i + 2047 <= #t) and (i + 2047) or #t
		res[#res + 1] = char(unpack(t, i, j))
	end
	return concat(res)
end

local function hex(s)
	local t = {}
	for i = 1, #s do t[i] = string.format("%02x", byte(s, i)) end
	return concat(t)
end

local function line(label, s)
	print(label .. " " .. hex(s))
end

local function numline(label, ...)
	print(label .. " " .. concat({ ... }, " "))
end

local function payload_seed(n, seed)
	local t = {}
	local x = seed
	for i = 1, n do
		x = band(x * 1664525 + 1013904223, 0xFFFFFFFF)
		t[i] = rshift(x, 24) % 256
	end
	return bytes_to_string(t)
end

local function xor_byte(s, pos, v)
	return sub(s, 1, pos - 1) .. char(bxor(byte(s, pos), v)) .. sub(s, pos + 1)
end

local function scenario_chromosomes()
	local data = payload_seed(200, 0xC0FFEE01)
	line("payload200", data)

	local c1 = chromosome.encode(7, data,
		{ gene_raw = 64, h = 3, units = 3, flags = 9 })
	if not c1 then print("chr_dense fail") return end
	line("chr_dense", c1)

	local c2 = chromosome.encode(7, data,
		{ gene_raw = 64, mode = "codon", h = 3, units = 3, flags = 9 })
	if not c2 then print("chr_codon fail") return end
	line("chr_codon", c2)

	local c3 = chromosome.set_generation(c1, 5)
	if not c3 then print("chr_gen5 fail") return end
	line("chr_gen5", c3)

	local c4 = chromosome.encode(7, data,
		{ gene_raw = 64, h = 3, units = 3, flags = 9, parity = 2 })
	if not c4 then print("chr_par2 fail") return end
	line("chr_par2", c4)

	local c5 = chromosome.encode(8, data,
		{ gene_raw = 48, mode = "codon", h = 3, units = 3, flags = 9, parity = 3 })
	if not c5 then print("chr_par3_codon fail") return end
	line("chr_par3_codon", c5)

	local r4 = chromosome.parse(c4)
	if not r4 then print("chr_par2_fields fail") return end
	numline("chr_par2_fields", r4.cen_version, r4.parity, r4.rawlen, r4.ngenes,
		#r4.genes)

	-- damage two data genes in the parity chromosome, then recover
	local broken = c4
	local rb = chromosome.parse(broken)
	for k = 1, 2 do
		broken = xor_byte(broken, rb.genes[k].offset + 12, 1)
	end
	local rb2 = chromosome.parse(broken)
	local damaged = 0
	for i = 1, #rb2.genes do
		if not rb2.genes[i].crc_ok then damaged = damaged + 1 end
	end
	numline("chr_par2_damaged", damaged, rb2.ngenes)
	local rec = chromosome.read(rb2)
	if rec then line("chr_par2_recover", rec) else print("chr_par2_recover fail") end

	-- the same corruption without parity cannot be recovered
	local plain = c1
	local rp = chromosome.parse(plain)
	for k = 1, 2 do
		plain = xor_byte(plain, rp.genes[k].offset + 12, 1)
	end
	local rp2 = chromosome.parse(plain)
	local rec2 = chromosome.read(rp2)
	print("chr_plain_recover " .. (rec2 and "ok" or "fail"))

	-- generation rewrite keeps the CEN2 centromere
	local c6 = chromosome.set_generation(c4, 7)
	if not c6 then print("chr_par2_gen7 fail") return end
	line("chr_par2_gen7", c6)

	-- parity-only reconstruction: erase every data gene, keep the parity
	local c7 = chromosome.encode(9, data,
		{ gene_raw = 64, h = 3, units = 3, flags = 9, parity = 4 })
	if not c7 then print("chr_par4 fail") return end
	local broken2 = c7
	local rb3 = chromosome.parse(broken2)
	for k = 1, 4 do
		broken2 = xor_byte(broken2, rb3.genes[k].offset + 12, 1)
	end
	local rb4 = chromosome.parse(broken2)
	local alive_par = 0
	for i = 1, #rb4.genes do
		local g = rb4.genes[i]
		if g.crc_ok and band(g.type, 2) == 2 then alive_par = alive_par + 1 end
	end
	numline("chr_par4_alive", alive_par, rb4.ngenes)
	local rec3 = chromosome.read(rb4)
	if rec3 then line("chr_par4_recover", rec3) else print("chr_par4_recover fail") end

	-- VIV14NB4NSH33: on-strand primer sites
	local c8 = chromosome.encode(11, data,
		{ gene_raw = 64, h = 3, units = 3, flags = 9, parity = 2, primer = 1234 })
	if not c8 then print("chr_banshee fail") return end
	line("chr_banshee", c8)
	local r8 = chromosome.parse(c8)
	if not r8 then print("chr_banshee_fields fail") return end
	numline("chr_banshee_fields", r8.cen_version, r8.parity, r8.primer,
		r8.primer_ok and 1 or 0, r8.ngenes, #r8.genes)

	-- destroy the reverse site: data stays readable, access is lost
	local pb = xor_byte(c8, #c8 - 23, 0x40)   -- first marker byte of the reverse site
	local r9 = chromosome.parse(pb)
	numline("chr_banshee_dark", r9.primer, r9.primer_ok and 1 or 0,
		chromosome.amplifiable(r9) and 1 or 0)
	local rr = chromosome.read(r9)
	if rr then line("chr_banshee_dark_read", rr) else print("chr_banshee_dark_read fail") end

	local c9 = chromosome.set_generation(c8, 11)
	if not c9 then print("chr_banshee_gen11 fail") return end
	line("chr_banshee_gen11", c9)
end

local function scenario_organisms()
	local d0 = payload_seed(300, 1)
	local d1 = payload_seed(180, 2)
	local e0 = payload_seed(300, 3)
	local e1 = payload_seed(180, 4)
	local o0 = { gene_raw = 64, h = 3, units = 3, flags = 0 }
	local o1 = { gene_raw = 64, h = 3, units = 3, flags = 5, parity = 2 }

	local org = organism.new({
		{ id = 0, data = d0, opts = o0 },
		{ id = 1, data = d1, opts = o1 },
	}, { max_gen = 40 })
	if not org then print("org fail") return end

	line("org1", organism.serialize(org))
	local s14 = organism.serialize14(org)
	local s14n = organism.serialize14n(org)
	line("org14", s14)
	line("org14n", s14n)
	numline("versions",
		organism.container_version(s14),
		organism.container_version(s14n),
		organism.container_version("VIV14NB4NSH33xx"))

	-- VIV14N roundtrip keeps the parity option
	local loaded = organism.deserialize(s14n)
	if loaded then
		line("org14n_rt", organism.serialize14n(loaded))
		numline("org14n_parity", loaded.chr_opts[2].parity)
	else
		print("org14n_rt fail")
	end

	-- VIV14 container carrying CEN2 strands infers the parity count
	local inf = organism.deserialize(s14)
	if inf then
		numline("org14_infer_parity", inf.chr_opts[2].parity)
	else
		print("org14_infer_parity fail")
	end

	-- cross before the damage/replication scenarios mutate org
	local orgB = organism.new({
		{ id = 0, data = e0, opts = o0 },
		{ id = 1, data = e1, opts = o1 },
	}, { max_gen = 40 })
	local child = organism.cross(org, orgB, 3)
	if child then
		line("org14n_cross", organism.serialize14n(child))
	else
		print("org14n_cross fail")
	end

	organism.damage(org, 3, 7)
	line("org14n_damage", organism.serialize14n(org))

	organism.replicate(org)
	line("org14n_replicate", organism.serialize14n(org))

	local m = organism.deserialize(s14n)
	if m then
		local mr = organism.mutate(m, 2, 5)
		if mr then
			numline("org_mutate_chr", mr.chr)
			line("org14n_mutate", organism.serialize14n(m))
		else
			print("org_mutate fail")
		end
	else
		print("org_mutate_load fail")
	end

	-- VIV14NB4NSH33 container
	local bopts = { o0, { gene_raw = 64, h = 3, units = 3, flags = 5, parity = 2, primer = 1234 } }
	local bo = organism.new({
		{ id = 0, data = d0, opts = bopts[1] },
		{ id = 1, data = d1, opts = bopts[2] },
	}, { max_gen = 40 })
	if bo then
		local bser = organism.serialize14nb(bo)
		if bser then
			line("org14nb", bser)
			local bl = organism.deserialize(bser)
			if bl then
				line("org14nb_rt", organism.serialize14nb(bl))
				numline("org14nb_opts", bl.chr_opts[2].parity, bl.chr_opts[2].primer)
			else
				print("org14nb_load fail")
			end
		else
			print("org14nb fail")
		end
	else
		print("org14nb_org fail")
	end
end

scenario_chromosomes()
scenario_organisms()
