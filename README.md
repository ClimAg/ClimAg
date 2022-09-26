# ClimAg

This research was funded by the [Environmental Protection Agency (EPA), Ireland][EPA]
project ["ClimAg: Multifactorial causes of fodder crises in Ireland and risks due to climate change"][ClimAg]
under the Climate Change Research Programme grant number 2018-CCRP-MS.50.

![ClimAg project logos](https://raw.githubusercontent.com/ClimAg/.github/main/images/logos.png)

## Features

WIP

## Licence

Copyright 2022 N. M. Streethran

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

  <https://www.apache.org/licenses/LICENSE-2.0>

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

## Credits

Pasture model ([ModVege]): This is the Python implementation of the ModVege
pasture model, translated from Java to Python by Chemin (2022).
The Java model was provided by R. Martin of INRAE UREP Clermont-Ferrand
for the original Python implementation.
The original ModVege pasture model was developed by
[Jouven et al. (2006)][Jouven].

Grass growth data for Northern Ireland has been extracted from bulletins of
the [GrassCheck] project by Agrisearch.

Some parts of this repository were created with [Cookiecutter] and the
[`audreyr/cookiecutter-pypackage`][audreyr] project template.

## References

### Software

- Chemin, Y. (2022). 'modvege', Python. [Online]. Available at
  <https://github.com/YannChemin/modvege> (Accessed 6 September 2022).

### Publications

- Jouven, M., Carrère, P. and Baumont, R. (2006). 'Model predicting dynamics
  of biomass, structure and digestibility of herbage in managed permanent
  pastures. 1. Model description', *Grass and Forage Science*, vol. 61, no. 2,
  pp. 112-124. DOI: [10.1111/j.1365-2494.2006.00515.x][Jouven].
- Huson, K. M., Lively, F. O., Aubry, A., Takahashi, T., Gordon, A. and McDonnell, D. A. (2020).
  'GrassCheck: monitoring grass growth and maximizing grass utilisation on UK farms',
  in Virkajärvi, P. et al. (eds), *Meeting the future demands for grassland production*,
  Grassland Science in Europe, Helsinki, Finland, European Grassland Federation,
  vol. 25, pp. 716–718. [Online]. Available at
  <https://www.europeangrassland.org/fileadmin/documents/Infos/Printed_Matter/Proceedings/EGF2020.pdf>
  (Accessed 13 September 2022).
- Hanrahan, L., Geoghegan, A., O'Donovan, M., Griffith, V., Ruelle, E.,
  Wallace, M. and Shalloo, L. (2017). 'PastureBase Ireland: A grassland
  decision support system and national database',
  *Computers and Electronics in Agriculture*, vol. 136, pp. 193–201.
  DOI: [10.1016/j.compag.2017.01.029][Hanrahan].

[EPA]: https://www.epa.ie/
[ClimAg]: https://www.ucc.ie/en/eel/projects/climag/
[ModVege]: https://github.com/YannChemin/modvege
[Jouven]: https://doi.org/10.1111/j.1365-2494.2006.00515.x
[Hanrahan]: https://doi.org/10.1016/j.compag.2017.01.029
[GrassCheck]: https://agrisearch.org/grasscheck
[Cookiecutter]: https://github.com/audreyr/cookiecutter
[audreyr]: https://github.com/audreyr/cookiecutter-pypackage
