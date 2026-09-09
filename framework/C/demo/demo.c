/* demo.c -- Interactive Bio-Tamagotchi: raise your data as an organism (C23) */
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vivi/organism.h"

// ANSI-цвета для красивого терминала
#define CLR_RESET   "\033[0m"
#define CLR_RED     "\033[1;31m"
#define CLR_GREEN   "\033[1;32m"
#define CLR_YELLOW  "\033[1;33m"
#define CLR_BLUE    "\033[1;34m"
#define CLR_CYAN    "\033[1;36m"
#define CLR_MAGENTA "\033[1;35m"

static void clear_screen(void)
{
#ifdef _WIN32
    // ANSI эскейп-код для очистки экрана в современных терминалах
    printf("\033[H\033[J");
#else
    printf("\033[H\033[2J");
#endif
}

static uint32_t get_seed(void)
{
    static uint32_t s = 1337;
    s ^= s << 13;
    s ^= s >> 17;
    s ^= s << 5;
    return s;
}

// Пул геномов для скрещивания с дикими особями
static const char *WILD_BRAINS[] = {
    "BRAIN:CURIOUS;SPEED:09;COLOR:GOLD;",
    "BRAIN:PREDATOR;SPEED:16;COLOR:PURPLE;",
    "BRAIN:PARASITE;SPEED:04;COLOR:BROWN;",
    "BRAIN:SWIFT;SPEED:22;COLOR:YELLOW;"
};
static const char *WILD_METABS[] = {
    "METABOLISM:PHOTOSYNTH;STAMINA:090;",
    "METABOLISM:SCAVENGER;STAMINA:140;",
    "METABOLISM:ANAEROBIC;STAMINA:080;",
    "METABOLISM:HYPERBOOST;STAMINA:220;"
};

static vivi_organism *create_wild_partner(const char **err)
{
    int idx = (int)(get_seed() % 4);
    const int ids[2] = { 0, 1 };
    const chr_opts co = { 1024, 0, 3, 4, 0 };
    const chr_opts cos[2] = { co, co };
    const uint8_t *data[2] = {
        (const uint8_t *)WILD_BRAINS[idx],
        (const uint8_t *)WILD_METABS[idx]
    };
    const size_t lens[2] = {
        strlen(WILD_BRAINS[idx]),
        strlen(WILD_METABS[idx])
    };

    vivi_organism *partner = nullptr;
    if (!organism_new(&partner, ids, cos, data, lens, 2, 60, err)) return nullptr;
    return partner;
}

int main(void)
{
    const char *err = nullptr;

    // Базовые опции хромосом
    const int ids[2] = { 0, 1 };
    const chr_opts co = { 1024, 0, 3, 4, 0 };
    const chr_opts cos[2] = { co, co };

    // Исходные данные (Адам)
    static const char AD0[] = "BRAIN:AGGRESSIVE;SPEED:12;COLOR:RED;";
    static const char AD1[] = "METABOLISM:CARNIVORE;STAMINA:100;";
    const uint8_t *ad[2] = { (const uint8_t *)AD0, (const uint8_t *)AD1 };
    const size_t dl[2] = { sizeof(AD0) - 1, sizeof(AD1) - 1 };

    // Эталонная стволовая клетка
    static const char ST0[] = "BRAIN:PRISTINE;SPEED:10;COLOR:WHITE;";
    static const char ST1[] = "METABOLISM:STEM_NICHE;STAMINA:150;";
    const uint8_t *st[2] = { (const uint8_t *)ST0, (const uint8_t *)ST1 };
    const size_t st_len[2] = { sizeof(ST0) - 1, sizeof(ST1) - 1 };

    vivi_organism *pet = nullptr;
    vivi_organism *stem_niche = nullptr;

    if (!organism_new(&pet, ids, cos, ad, dl, 2, 60, &err) ||
        !organism_stem(&stem_niche, ids, cos, st, st_len, 2, 60, &err)) {
        fprintf(stderr, "Initialization failed: %s\n", err);
        return 1;
    }
    organism_attach_stem(pet, stem_niche);

    char last_event[256] = "Organism born into the primordial soup.";

    while (true) {
        clear_screen();

        // 1. Считываем состояние организма
        vivi_bytes *karyotype = nullptr;
        cell_report rep = { 0 };
        bool is_alive = organism_read(&karyotype, pet, &rep, &err);

        // 2. Отрисовка UI
        printf(CLR_CYAN "========================================================================\n" CLR_RESET);
        printf(CLR_CYAN "            VIVIAN - SYNTHETIC LIFE SIMULATOR (C23 ENGINE)             \n" CLR_RESET);
        printf(CLR_CYAN "========================================================================\n" CLR_RESET);

        // Лицо и статус существа
        if (!is_alive) {
            printf("\n   " CLR_RED "( X _ X )   STATUS: COLLAPSED / DEAD" CLR_RESET "\n");
        } else if (rep.renewed) {
            printf("\n   " CLR_MAGENTA "( * v * )!  STATUS: REGENERATED FROM STEM NICHE" CLR_RESET "\n");
        } else if (pet->generation >= pet->max_gen) {
            printf("\n   " CLR_YELLOW "( - _ - )   STATUS: SENESCENT (HAYFLICK LIMIT REACHED)" CLR_RESET "\n");
        } else {
            printf("\n   " CLR_GREEN "( ^ _ ^ )b  STATUS: HEALTHY & ACTIVE" CLR_RESET "\n");
        }

        // Шкала поколений (Лимит Хейфлика)
        printf("\n  " CLR_YELLOW "Age / Generation: " CLR_RESET "[%d / %d] ", pet->generation, pet->max_gen);
        printf("[");
        for (int i = 0; i < 30; ++i) {
            if (i < (pet->generation * 30 / pet->max_gen)) printf(CLR_RED "|" CLR_RESET);
            else printf(CLR_GREEN "." CLR_RESET);
        }
        printf("]  Stem Link: %s\n", pet->stem_source ? CLR_GREEN "ONLINE" CLR_RESET : CLR_RED "OFFLINE" CLR_RESET);

        // Отображение хромосом
        printf("\n" CLR_CYAN "--- KARYOTYPE & EXPRESSED TRAITS ---------------------------------------" CLR_RESET "\n");
        if (is_alive && karyotype) {
            for (size_t i = 0; i < pet->nchr; ++i) {
                printf("  Chr %d: " CLR_GREEN "%.*s" CLR_RESET "\n",
                    pet->chr_ids[i], (int)karyotype[i].len, (const char *)karyotype[i].data);
            }
            vivi_bytes_free_n(karyotype, pet->nchr);
        } else {
            printf(CLR_RED "  [GENOME DESTROYED - HOMOLOGOUS REPAIR FAILED]\n" CLR_RESET);
        }

        // Лог последнего события
        printf("\n" CLR_CYAN "--- EVENT LOG ----------------------------------------------------------" CLR_RESET "\n");
        printf("  >> %s\n", last_event);

        // Меню действий
        printf("\n" CLR_CYAN "--- ACTIONS ------------------------------------------------------------" CLR_RESET "\n");
        printf("  [1] " CLR_YELLOW "Zap with Cosmic Rays" CLR_RESET "   (Test Diploid Homologous Repair)\n");
        printf("  [2] " CLR_MAGENTA "Lethal Gamma Burst"  CLR_RESET "   (Blast with 30 mutations -> Stem Rescue)\n");
        printf("  [3] " CLR_GREEN "Induce Point Mutation" CLR_RESET " (Evolutionary single-bit flip with Chaskey)\n");
        printf("  [4] " CLR_CYAN "Sex / Crossing-Over"   CLR_RESET "   (Breed with wild mate -> new child)\n");
        printf("  [5] " CLR_YELLOW "Mitosis / Cellular Age" CLR_RESET " (Divide cell -> advance Hayflick +1)\n");
        printf("  [6] " CLR_GREEN "Stem Rejuvenation"     CLR_RESET "   (Reset generation back to 0)\n");
        printf("  [7] " CLR_BLUE "Save to pet.viv1"      CLR_RESET "    (Binary serialization)\n");
        printf("  [8] " CLR_BLUE "Load from pet.viv1"    CLR_RESET "  (Binary deserialization)\n");
        printf("  [0] Exit Simulator\n");
        printf("\nSelect action > ");

        int choice = -1;
        if (scanf("%d", &choice) != 1) {
            int c;
            while ((c = getchar()) != '\n' && c != EOF) {}
            continue;
        }

        if (choice == 0) break;

        switch (choice) {
            case 1: { // Умеренная радиация
                (void)organism_damage(pet, 2, get_seed(), &err);
                vivi_bytes *dummy = nullptr;
                cell_report r;
                (void)organism_read(&dummy, pet, &r, &err);
                if (dummy) vivi_bytes_free_n(dummy, pet->nchr);
                snprintf(last_event, sizeof(last_event),
                    "Radiation (2 hits): Repaired: %d, Structural: %d, Dead: %d",
                    r.repaired, r.structural, r.dead);
                break;
            }
            case 2: { // Смертельная радиация
                (void)organism_damage(pet, 30, get_seed(), &err);
                vivi_bytes *dummy = nullptr;
                cell_report r;
                bool ok = organism_read(&dummy, pet, &r, &err);
                if (dummy) vivi_bytes_free_n(dummy, pet->nchr);
                if (ok && r.renewed) {
                    snprintf(last_event, sizeof(last_event),
                        CLR_MAGENTA "LETHAL DOSE (30 hits)! Internal repair failed, but STEM NICHE restored the organism!" CLR_RESET);
                } else {
                    snprintf(last_event, sizeof(last_event),
                        CLR_RED "LETHAL DOSE! Organism perished permanently." CLR_RESET);
                }
                break;
            }
            case 3: { // Точечная мутация
                org_mut_result mu;
                if (organism_mutate(&mu, pet, 1, get_seed(), &err)) {
                    snprintf(last_event, sizeof(last_event),
                        "Viable mutation accepted on Chr %d! Chaskey MAC updated.", mu.chr);
                    vivi_bytes_free(&mu.data);
                } else {
                    snprintf(last_event, sizeof(last_event), "Mutation rejected: %s", err ? err : "unknown");
                }
                break;
            }
            case 4: { // Скрещивание (Мейоз / Кроссинговер)
                vivi_organism *mate = create_wild_partner(&err);
                if (mate) {
                    vivi_organism *child = nullptr;
                    if (organism_cross(&child, pet, mate, get_seed(), &err)) {
                        organism_attach_stem(child, stem_niche);
                        organism_free(pet);
                        pet = child;
                        snprintf(last_event, sizeof(last_event),
                            CLR_GREEN "Bred with wild mate! Child born with mosaic genome. Generation: 0." CLR_RESET);
                    } else {
                        snprintf(last_event, sizeof(last_event), "Breeding failed: %s", err ? err : "unknown");
                    }
                    organism_free(mate);
                } else {
                    snprintf(last_event, sizeof(last_event), "Wild mate creation failed: %s", err ? err : "unknown");
                }
                break;
            }
            case 5: { // Митоз (деление и старение)
                vivi_organism *daughter = nullptr;
                if (organism_mitosis(&daughter, pet, &err)) {
                    organism_free(pet);
                    pet = daughter;
                    snprintf(last_event, sizeof(last_event),
                        "Mitosis successful! Organism divided. Generation advanced to %d.", pet->generation);
                } else {
                    snprintf(last_event, sizeof(last_event),
                        CLR_RED "Mitosis blocked: %s" CLR_RESET, err ? err : "senescent");
                }
                break;
            }
            case 6: { // Омоложение через стволовые клетки
                if (organism_renew(pet, stem_niche, &err)) {
                    snprintf(last_event, sizeof(last_event),
                        CLR_CYAN "Rejuvenation complete! Generation counter reset to 0." CLR_RESET);
                } else {
                    snprintf(last_event, sizeof(last_event), "Rejuvenation failed: %s", err ? err : "unknown");
                }
                break;
            }
            case 7: { // Сохранение в бинарник VIV1
                vivi_bytes saved;
                if (organism_serialize(&saved, pet, &err)) {
                    FILE *f = fopen("pet.viv1", "wb");
                    if (f) {
                        fwrite(saved.data, 1, saved.len, f);
                        fclose(f);
                        snprintf(last_event, sizeof(last_event),
                            "Saved genome to pet.viv1 (%zu bytes, ~3%% overhead).", saved.len);
                    } else {
                        snprintf(last_event, sizeof(last_event), "File write error.");
                    }
                    vivi_bytes_free(&saved);
                }
                break;
            }
            case 8: { // Загрузка из pet.viv1
                FILE *f = fopen("pet.viv1", "rb");
                if (f) {
                    fseek(f, 0, SEEK_END);
                    long sz = ftell(f);
                    fseek(f, 0, SEEK_SET);
                    uint8_t *buf = malloc(sz);
                    if (buf && fread(buf, 1, sz, f) == (size_t)sz) {
                        vivi_organism *loaded = nullptr;
                        if (organism_deserialize(&loaded, buf, sz, &err)) {
                            organism_free(pet);
                            pet = loaded;
                            organism_attach_stem(pet, stem_niche);
                            snprintf(last_event, sizeof(last_event),
                                "Loaded genome from pet.viv1 successfully!");
                        } else {
                            snprintf(last_event, sizeof(last_event), "Corrupt container: %s", err);
                        }
                    }
                    free(buf);
                    fclose(f);
                } else {
                    snprintf(last_event, sizeof(last_event), "File pet.viv1 not found. Save first!");
                }
                break;
            }
            default:
                snprintf(last_event, sizeof(last_event), "Unknown command.");
                break;
        }
    }

    organism_free(pet);
    organism_free(stem_niche);
    printf("\nSimulator terminated. All biological structures safely deallocated.\n");
    return 0;
}



