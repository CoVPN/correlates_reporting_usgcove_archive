source(here::here("code", "params.R"))

# 1. Baseline risk factors
# 2. Baseline risk factors and the bAb anti-Spike markers,  
# 3. Baseline risk factors and the bAb anti-RBD markers,  
# 4. Baseline risk factors and the pseudovirus-nAb ID50 markers,   
# 5. Baseline risk factors and the pseudovirus-nAb ID80 markers,   
# 6. Baseline risk factors and the live virus-nAb MN50 markers,   
# 7. Baseline risk factors and the bAb markers and the pseudovirus-nAb ID50markers,   
# 8. Baseline risk factors and the bAb markers and the pseudovirus-nAb ID80markers,   
# 9. Baseline risk factors and the bAb markers and the live virus-nAb MN50markers,    
# 10. Baseline risk factors and the bAb markers and the combination scores across the five markers [PCA1, PCA2, FSDAM1/FSDAM2 (the first two components of nonlinear PCA), and the maximum signal diversity score (He and Fong (2019)],   
# 11. Baseline risk factors and all individual marker variables,    
# 12. Baseline risk factors and all individual marker variables and all combination scores (full model)"

bindSpike_markers <- get_markers("bindSpike")
bindRBD_markers <- get_markers("bindRBD")
pseudoneutid50_markers <- get_markers("pseudoneutid50")
pseudoneutid80_markers <- get_markers("pseudoneutid80")
liveneutmn50_markers <- get_markers("liveneutmn50")

set1 <- bl_demo_var
set2 <- c(bl_demo_var, bindSpike_markers)
set3 <- c(bl_demo_var, bindRBD_markers)
set4 <- c(bl_demo_var, pseudoneutid50_markers)
set5 <- c(bl_demo_var, pseudoneutid80_markers)
set6 <- c(bl_demo_var, liveneutmn50_markers)
set7 <- c(bl_demo_var, bindSpike_markers, bindRBD_markers, pseudoneutid50_markers)
set8 <- c(bl_demo_var, bindSpike_markers, bindRBD_markers, pseudoneutid80_markers)
set9 <- c(bl_demo_var, bindSpike_markers, bindRBD_markers, liveneutmn50_markers)
set10 <- c(bl_demo_var, bindSpike_markers, bindRBD_markers, "pc1", "pc2", "max_sig")
set11 <- c(bl_demo_var, bindSpike_markers, bindRBD_markers, 
           pseudoneutid50_markers, pseudoneutid80_markers, liveneutmn50_markers)
set12 <- c(bl_demo_var, bindSpike_markers, bindRBD_markers, 
           pseudoneutid50_markers, pseudoneutid80_markers, liveneutmn50_markers,
           "pc1", "pc2", "max_sig")

full_model <- set12