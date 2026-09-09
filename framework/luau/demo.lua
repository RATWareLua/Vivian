--!native
local organism = require("./organism")

print("=== 1. CREATING PARENTS & STEM CELL ===")
local adam = organism.new({
    { id = 0, data = "BRAIN:AGGRESSIVE;SPEED:12;COLOR:RED;", opts = { mode = "dense" } },
    { id = 1, data = "METABOLISM:CARNIVORE;STAMINA:100;",  opts = { mode = "dense" } },
})

local eve = organism.new({
    { id = 0, data = "BRAIN:DEFENSIVE;SPEED:08;COLOR:BLUE;", opts = { mode = "dense" } },
    { id = 1, data = "METABOLISM:HERBIVORE;STAMINA:150;",  opts = { mode = "dense" } },
})

-- Создаем эталонную стволовую клетку для регенерации
local stem_cell = organism.stem({
    { id = 0, data = "BRAIN:STANDARD;SPEED:10;COLOR:GREEN;", opts = { mode = "dense" } },
    { id = 1, data = "METABOLISM:OMNIVORE;STAMINA:120;",    opts = { mode = "dense" } },
})

print("=== 2. SEXUAL REPRODUCTION (CROSSING OVER) ===")
local child = organism.cross(adam, eve, 42)
organism.attach_stem(child, stem_cell) -- прикрепляем стволовую нишу
print("Child born! Generation:", child.generation)

print("\n=== 3. HOMOLOGOUS REPAIR (1 hit per chr, loci differ) ===")
-- Наносим легкий урон: гомологи повреждаются в разных местах
organism.damage(child, 1, 101)

local data, rep_or_err, rep = organism.read(child)
local r = rep or rep_or_err
print(string.format("Checkpoint: Repaired: %d, Structural: %d, Dead: %d",
    r.repaired, r.structural, r.dead))
print("Chr 0:", data[0])
print("Chr 1:", data[1])

print("\n=== 4. VIABLE MUTATION (EVOLUTIONARY LEAP) ===")
local mut_res = organism.mutate(child, 1, 555)
local mut_data = organism.read(child)
print(string.format("Mutated Chr %d successfully! New data:", mut_res.chr))
print(" ->", mut_data[mut_res.chr])

print("\n=== 5. SERIALIZATION (VIV1 CONTAINER) ===")
local saved = organism.serialize(child)
print(string.format("Serialized container size: %d bytes", #saved))

local loaded = organism.deserialize(saved)
local loaded_data = organism.read(loaded)
print("Loaded from disk Chr 1:", loaded_data[1])

print("\n=== 6. LETHAL RADIATION & STEM CELL RESCUE ===")
print("Blasting with 30 mutations (both homologs destroyed)...")
organism.damage(child, 30, 777)

-- При смертельной дозе восстановить по гомологу нельзя,
-- но клетка обращается к stem_source и возрождается!
local rescue_data, rescue_rep = organism.read(child)
if rescue_data then
    print("Organism regenerated from STEM NICHE! Renewed:", rescue_rep.renewed or false)
    print("Restored Chr 0:", rescue_data[0])
else
    print("Organism perished completely.")
end
