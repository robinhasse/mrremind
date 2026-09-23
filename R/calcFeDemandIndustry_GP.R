calcFeDemandIndustry_GP <- function(scenarios) {
  scenarioIndustry <- setdiff(scenarios, "SSP2_GP")

  feIndustry <- calcOutput("FeDemandIndustry",
                           scenarios = scenarioIndustry,
                           warnNA = FALSE,
                           aggregate = FALSE)


  # Copy SSP2 to good performance scenario for industry FE demand
  feIndustry <- mbind(
    feIndustry[, , "SSP2", invert = TRUE],
    setItems(feIndustry[, , "SSP2"], 3.1, "SSP2_GP")
  )

  return(list(x = feIndustry,
              description = "demand pathways for final energy demand in industry",
              unit = "EJ, except ue_cement (Gt), ue_primary_steel and ue_secondary_steel (Gt) and ue_chemicals and ue_otherInd ($tn)"))
}
