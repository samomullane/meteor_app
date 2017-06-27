# meteor_app:
## NASA JPL meteor and comet data visualization

The finished app can be found [here](https://samomullane.shinyapps.io/meteor_app/).
A full write up of this project can be accessed via [this link](http://blog.nycdatascience.com/student-works/r-shiny/meteor_impact/).

The contents of this main folder are the shiny R source files (global, server and ui files). 
The folders are img (images used in shiny application) and cleaned_data (data used for analysis in app and blog post).

The key packages used for this project are:
1. Shiny - overall framework
2. dplyr - data manipulation and preparation
3. plotly - prettifying ggplots and making them interactive
4. leaflet - showing the city-sized craters created by meteor impacts
