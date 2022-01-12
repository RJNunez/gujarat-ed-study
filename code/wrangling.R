# -- Set up
source("code/init.R")

# -- Data paths
paths <- sort(list.files("data/csvs", full.names = TRUE))
idx   <- 1:length(paths)

# -- Loading and wrangling
dat <- map_df(idx[-91], function(i){
  
  # -- Index for columns of interest
  if(i == 16) {
    col_idx <- c(1, 2, 5, 6, 7, 8)
  } else {
    col_idx <- c(2, 3, 6, 7, 8, 9)#c(1, 3, 6, 7, 8, 9)
  }

  # -- Wrangling data
  dat <- read_csv(paths[i]) %>%
    dplyr::select(col_idx) %>%
    setNames(c("subdistrict", "village", "date", "gender", "age", "place")) %>%
    fill(subdistrict, village) %>%
    mutate(i = i, 
           date    = as.Date(date, format = "%d/%m/%Y"),
           village = str_to_title(village),
           village = str_replace_all(village, "\\(M\\)", ""),
           village = str_replace_all(village, " ", ""),
           tmp2    = case_when(village == "Dehgam"   ~ "Dahegam",
                               village == "DwarkaNp" ~ "Dwarka",
                               village == "Visadar"  ~ "Visavadar",
                               village == "DhaneraNagarPalika" ~ "Dhanera",
                               village == "TarsadiNagarpalika" ~ "Tarsadi",
                               village %in% c("Savrakundla", "Savrakundala") ~ "Savarkundla",
                               village %in% c("Vallabhaveedhanagar", "Vallabhaveedyanagar", "VallabhVeedyanagar", "VallbhVeedhanagar") ~ "Vallabh Vidhyanagar",
                               village %in% c("JadanNagarpalika", "NagarpalikaJasdan") ~ "Jasdan",
                               village %in% c("Navsari_np", "Navsari-VijalporNagarpalika") ~ "Navsari",
                               village == "Gadhda"      ~ "Gadhada",
                               village == "Khedbrhma"   ~ "Khedbrahma",
                               village == "KillaPardi"  ~ "Pardi",
                               village == "Kuteeyana"   ~ "Kutiyana",
                               village == "Vakaner"     ~ "Wankaner",
                               village == "Lunavada"    ~ "Lunawada",
                               village == "Dhangdhra"   ~ "Dhrangadhra",
                               village == "Kaalol"      ~ "Kalol",
                               village == "Khabhalia"   ~ "Khambhalia",
                               village == "Mandavi"     ~ "Mandvi",
                               village == "Barvala"     ~ "Barwala",
                               village == "Rajpipala"   ~ "Rajpipla",
                               village == "Savali"      ~ "Savli",
                               village == "Disa"        ~ "Deesa",
                               village == "Jafarabad"   ~ "Jafrabad",
                               TRUE ~ village),
           subdistrict = str_to_title(subdistrict),
           tmp      = case_when(subdistrict == "Arvali" ~ "Arvalli",
                                subdistrict == "Khabhalia"  ~ "Khambhalia",
                                subdistrict == "Banaskntha" ~ "Banaskantha",
                                subdistrict == "Barvala" ~ "Barwala",
                                subdistrict %in% c("Devbumi Dwarka", "Devbhoomi Dwarka") ~ "Devbhumi Dwarka",
                                subdistrict == "Dhangdhra" ~ "Dhrangadhra",
                                subdistrict == "Savali" ~ "Savli",
                                TRUE ~ subdistrict),
           gender   = case_when(gender %in% c("Male", "M", "m", "M .", "MAIL", "M #", "Ml", "Me", "M ■", "Mae", "M »", "M *", "M j", "«m", "M '") ~ "male",
                                gender %in% c("Female", "F", "f", "F ~1", "FF", "(F", "FEMAIL", "F t", "F ‘", "F.", "FM", "F •", "F .") ~ "female",
                                TRUE ~ "other"),
           place = case_when(str_detect(place, "Medicare|SPECIALITY|speciality|CENTOR|OFFICE|I.C.U|Lifecare|Nursing|NURSING|MEDICA|Center|center|Heart|DZ|COVID|CENTER|CENTRE|Centre|m^dicai|Cardiolog|Dr|Doctor|MEDICAL|C .H .C .|C H C|C H.C|C.H.C|C.H..C.|DOCTOR|ICU|narsing|nursing|DR|medical|centre|CHC|H.O|CARE|Hoa|Hoe|Hot|hosp|HOSP|Hosp|al") ~ "Medical Center",
                             str_detect(place, "House|Home|Hous|Hoitse|H0U86|H0U88|Heuser|Hgusg|Hotise|house|Hou’se|Hou6e|HOU86|HOU88|HoU8B|Hcuse|Iouse|-Iouse|1 iouse|HOME|HOUSE|Hoose|-iouse|4ouse|-louse|louse|Hours|Houee|Mouse|ilouse|ouse") ~ "House",
                             str_detect(place, "Other") ~ "Other"),
           # place2 = case_when(str_detect(place, "ther") ~ "Other",
           #                    str_detect(place, "Medicare|SPECIALITY|speciality|CENTOR|OFFICE|I.C.U|Lifecare|Nursing|NURSING|MEDICA|Center|center|Heart|DZ|COVID|CENTER|CENTRE|Centre|m^dicai|Cardiolog|Dr|Doctor|MEDICAL|C .H .C .|C H C|C H.C|C.H.C|C.H..C.|DOCTOR|ICU|narsing|nursing|DR|medical|centre|CHC|H.O|CARE|Hoa|Hoe|Hot|hosp|HOSP|Hosp|al") ~ "Medical Center",
           #                    TRUE ~ "House"),
           age      = str_replace(age, "Year", ""),
           age      = as.numeric(age),
           age      = cut(age, breaks = seq(0, 100, by = 10), ordered_result = TRUE)
    ) %>%
    dplyr::select(-subdistrict, -village) %>%
    rename(subdistrict = tmp,
           village  = tmp2)
})

# -- Saving data
save(dat, file = "data/rdas/municipality-data.rda", compress = "xz")

