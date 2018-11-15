setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

#Reading and tidying----
source('2.Libraries.R')
source('3.Functions.R')
source('4.TrmmLoading.R')
source('5.LandfluxAetLoading.R')
source('6.PrincetonPetLoading.R')
source('7.LoadCru.R')
source('8.ComputeSpei.R')
source('9.WaterDeficitTibbleCompilation.R')
source('10.CwdAnomalies.R')

#Plotting
source('11.DryMonthsAndAnnualRainfallPlotting.R')
source('12.GatheringAndPlottingAnomalies.R')
source('13.Correlation-GLK.R')
source('14.TemporalEvolutionOfDroughtExtent.R')
source('15.SpatialCorrelationAcrossIndices.R')
source('16.spiSpei2005.R')
source('17.Annexes.R')