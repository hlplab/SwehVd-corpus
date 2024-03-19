# Constants ------------------------------------------------------------------
theme_set(
  theme_bw(base_size = 11) + 
    theme(
      panel.grid = element_blank()))

levels.condition.wade = c("native mean & correlation, minimal variance", "native mean & correlation, native variance", "native mean & correlation, nonnative variance", "nonnative mean & correlation, minimal variance", "nonnative mean & correlation, native variance", "nonnative mean & correlation, nonnative variance")
levels.vowel.IPA = c("[i]", "[ɪ]", "[ɛ]", "[æ]", "[ʌ]", "[ʊ]", "[u]", "[ɑ]")
levels.vowel.Arpabet = c("iy1", "ih1", "eh1", "ae1", "ah1", "uh1", "uw1", "aa1")
levels.vowel.SweFA = c("ii1", "ih1", "yy1", "yh1", "uu1", "uh1", "ee1", "eh1", "ae1", "{:", "{", "oe1", "oeh1", "9:", "9", "aa1", "ah1", "oa1", "oah1", "oo1", "oh1")
levels.vowel.IPA.swe = c("[iː]", "[ɪ]", "[yː]", "[ʏ]", "[ʉː]", "[ɵ]", "[eː]", "[ɛ]", "[ɛː]", "[æː]", "[æ]", "[øː]", "[ø]", "[œː]", "[œ]", "[ɑː]", "[a]", "[oː]", "[ɔ]", "[uː]", "[ʊ]")
levels.vowel.IPA.swe.long = c("[iː]", "[yː]", "[ʉː]", "[eː]", "[ɛː]", "[æː]", "[øː]", "[œː]", "[ɑː]", "[oː]", "[uː]")
levels.vowel.IPA.swe.long.noallo = c("[iː]", "[yː]", "[ʉː]", "[eː]", "[ɛː]", "[øː]", "[ɑː]", "[oː]", "[uː]")
levels.vowel.IPA.swe.short = c("[ɪ]", "[ʏ]", "[ɵ]", "[ɛ]", "[æ]", "[ø]", "[œ]", "[a]", "[ɔ]", "[ʊ]")
levels.vowel.SweFA.wrongEncod = c("ii1", "ih1", "yy1", "yh1", "uu1", "uh1", "ee1", "eh1", "ae1", "{:", "{", "oe1", "oeh1", "09:00", "9", "aa1", "ah1", "oa1", "oah1", "oo1", "oh1")

levels.response.natural.swe.long <- c("hid", "hyd", "hud", "hed", "häd", "härd", "höd", "hörd", "had", "håd", "hod")
levels.response.natural.swe.long.noallo <- c("hid", "hyd", "hud", "hed", "häd", "höd", "had", "håd", "hod")

# Color codes for plotting
colors.vowel <- c("#1B9E77", "#D95F02", "#7570B3", "#FB9A99", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#1F78B4", "#33A02C")
names(colors.vowel) <- levels.vowel.IPA
colors.vowel.word <- c("#1B9E77", "#D95F02", "#7570B3", "#FB9A99", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666", "#1F78B4", "#33A02C")
names(colors.vowel.word) <- levels.response.natural.swe.long
colScale <-   scale_colour_manual(name = "Vowel",values = colors.vowel)
colors.vowel.swe <- c("#66C2A5", "#1B9E77", "#FC8D62", "#D95F02", "#8DA0CB", "#7570B3", "#FB9A99", "#E78AC3", "#E7298A", "#A6D854", "#66A61E", "#FFD92F", "#E6AB02", "#E5C494", "#A6761D", "#B3B3B3", "#666666", "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C")
names(colors.vowel.swe) <- levels.vowel.IPA.swe

vowelPlot_components <- list(
  axis = c(scale_x_reverse(), scale_y_reverse()),
  scale = colScale, 
  theme(legend.position="right"), 
  guides(color = guide_legend(order=1), shape = guide_legend(order=2)))
