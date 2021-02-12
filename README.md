# change-metabolomics

*Publication coming soon*

Growing industry and expanding populations of people have led to big changes in the global nitrogen cycle. For example, a lot of the nitrogen put down as fertilizer and produced from burning fossil fuels ends up in natural environments. I wanted to know how a natural environment I'm passionate about, the shortgrass steppe of Western North America, responds to increased nitrogen. Specifically, what's going on at the micro level? This means targeting cellular metabolites and understanding movement of water in individual plants (aka photosynthesis).

![Shortgrass steppe](https://upload.wikimedia.org/wikipedia/commons/thumb/0/00/Llano_Estacado_Caprock_Escarpment_south_of_Ralls_TX_2009.jpg/1920px-Llano_Estacado_Caprock_Escarpment_south_of_Ralls_TX_2009.jpg)  

*Shortgrass steppe*

I chose to look at two really common, or dominant, species in the shortgrass steppe to answer this question. Using metabolite ("metabolomics") data, I found that blue grama grass (*Bouteloua gracilis*) and scarlet globemallow (*Sphaeralcea coccinea*) differed a LOT at the metabolite level. When we added nitrogen, we found a much smaller effect. Approximately 8% of metabolites clustered into a network responded to nitrogen. I also used [maximal information coefficients](https://en.wikipedia.org/wiki/Maximal_information_coefficient) to determine which metabolites responded in a nonlinear way to nitrogen. None of these metabolites overlapped in the two species, indicating different responses to nitrogen at the metabolic level.. although in the big scheme of things this effect is small!

Interestingly, blue grama increased it's photosynthetic rate with more nitrogen, while globemallow did not. This *could* mean that in a future with more nitrogen, blue grama grows faster and becomes more abundant! However, there are other plants in the community to consider as competitors, such as fast growing wild rye (*Elymus elymoides*). We'll have to continue monitoring the shortgrass steppe responses to nitrogen for the clearest picture.

![Plant metabolites and PCs](https://github.com/avahoffman/change-metabolomics/blob/master/figures/PCA_MIC.png)  

*(a): Principal components analysis showing large separation between blue grama grass (*Bouteloua gracilis*) and scarlet globemallow (*Sphaeralcea coccinea*) metabolites. (b): maximal information coefficient analysis revealed several metabolites with a significant (p<0.001) positive relationship and (c) negative relationship with nitrogen. Each line represents an individual metabolite in (b) and (c). Hashed lines represent metabolites also responded to nitrogen more generally*

To run this analysis, first clone the repository. Navigate to the data (stored on Figshare) and ensure it's added to the data directory within your cloned repo. Then, you'll want to navigate to the run file:
```
├── src
│   ├── run.R
```
Open directly into RStudio (or other IDE) to ensure the `here` package works as expected.
