# Yukon-Chinook-CU-reconstructions
Code and data to reconstruct Yukon Chinook salmon Conservation Units over time based on genetics of fish sampled at the Alaska-Yukon Territory border, annual estimates of border passage and reconstructed harvest rates for the Canadian-origin stock aggregate. 

Details on the run-reconstruction and spawner-recruitment analyses are available in: 
> Connors, B.M., O’Dell, A., Hunter, H., Glaser, D., Gill, J., Rossi, S., and Churchland, C. 2025. Stock status and biological and fishery consequences of alternative harvest and rebuilding actions for Yukon River Chinook salmon (*Oncorhynchus tshawytscha*). DFO Can. Sci. Advis. Sec. Res. Doc. 2025/nnn. iv + 130 p.

and in particular Appendices A and B. The code and data to reproduce the above report are in this [Github repository](https://github.com/Pacific-salmon-assess/Yukon-Chinook-stock-assessment-2025/tree/main#yukon-chinook-stock-assessment-2025). 

### This repository is for annual updates to the run-reconstructions. 
To generate updated reconstructions:
 - Run the run-reconstruction model fitting code ([`RR.R`](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/run-reconstructions/RR.R)) to fit the multi-Conservation Unit border passage reconstructions model, and generate time series of spawners and harvest.
 - Run the spawner-recruit model fitting code ([`SR_fit.R`](https://github.com/Pacific-salmon-assess/yukon-CK-ResDoc/blob/main/analysis/R/SR_fit.R)) to fit the spawner-recruit models, and generate estimates of biological reference points.

