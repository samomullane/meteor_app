library(data.table)
library(magrittr)
library(ggplot2)
#library(globe)
#globeearth(eye=place("newyorkcity"))
install.packages('shinyRGL')
library(shinyRGL)
library(dplyr)


setwd('~/Desktop/nycdsa/shiny_comet/')

#############################################################
#read in sentry database for probable impacts
sentry_dt <- fread('./sentry.csv')

sentry_dt[, c('Palermo.Scale.min') := sentry_dt[,8]]
sentry_dt[, c('Palermo.Scale.max') := sentry_dt[,9]]
sentry_dt[,9] = NULL
sentry_dt[,8] = NULL
names(sentry_dt)

#Rename columns in sentry dt to space-less, parenthesis-free

names(sentry_dt) %<>% 
  gsub(pattern='\\s\\((.*?)\\)', 
       x=., 
       replacement = '\\2') %>% 
  gsub(pattern='\\s',
       x=.,
       replacement = '\\.\\1')

sentry_dt$Object.Designation.. <- gsub(pattern='(.?)*\\((.*?)\\)', 
     x=sentry_dt$Object.Designation.., 
     replacement = '\\2')

sentry_dt$Year.Min <- sentry_dt$Year.Range.. %>% 
  gsub(pattern='\\-[0-9]*',
     x=.,
     replacement = '\\1')

sentry_dt$Year.Max <- sentry_dt$Year.Range.. %>% 
  gsub(pattern='[0-9]*\\-',
     x=.,
     replacement = '\\1')

sentry_dt <- sentry_dt[, lapply(.SD, as.numeric), by=c('Object.Designation..',
                                          'V11',
                                          'Year.Range..')]

#############################################################
#read in small body database (main)

small_body_dt <- fread('./small_body.csv')

#Take out columns with almost entirely NAs
small_body_dt <- small_body_dt[, 
                             -c('G', 'A1', 'A2', 'A3', 'DT',
                                'sigma_e', 'sigma_a',
                                'sigma_q', 'sigma_i',
                                'sigma_om', 'sigma_w',
                                'sigma_ma', 'sigma_ad',
                                'sigma_n'),
                             with = F]

small_body_dt$full_name <- gsub(pattern='(.?)*\\((.?)*\\)',
     x=small_body_dt$full_name, 
     replacement = '\\2') %>% 
  gsub(pattern='\\s+\\s',
       x=., 
       replacement = '')

#Class level statistics
sbdt_summary <- small_body_dt %>% 
  group_by(class) %>% 
  summarise(class_count = n(),
            avg_H = mean(H, na.rm=T),
            avg_e = mean(e),
            avg_a = mean(a),
            avg_q = mean(q),
            avg_i = mean(i),
            avg_alb = mean(albedo),
            avg_diam = mean(diameter, na.rm=T),
            avg_moid = mean(moid)
  )

sbdt_summary$class = as.factor(sbdt_summary$class)

ggplot(sbdt_summary, aes(x = factor(1),
                         fill = factor(class_count))) +
  geom_bar(width=1) + coord_polar(theta = 'y')

small_body_dt$class <- as.factor(small_body_dt$class)

ggplot(small_body_dt, aes(x = diameter, y = H)) + geom_jitter() + facet_wrap(~class)

ggplot(sentry_dt, aes(x = Estimated.Diameter, y = H)) + geom_jitter(aes(color=Vinfinity))

#cumulative probability that at least one object impacts earth:
1-cumprod(1-sentry_dtc$Impact.Probability)[nrow(sentry_dt)]
# result = 8.694%

#Full data from both tables is available on 9 entries
small_body_join <- inner_join(sentry_dt, small_body_dt, by=c("Object.Designation.."="full_name"))

#Filter out all rows with NAs
small_body_dt <- small_body_dt[,-c('GM'), with=F]
small_body_dt <- small_body_dt[complete.cases(small_body_dt$a),]
small_body_dt <- small_body_dt[complete.cases(small_body_dt$per_y),]
small_body_dt <- small_body_dt[complete.cases(small_body_dt$ma),]
small_body_dt <- small_body_dt[complete.cases(small_body_dt$om),]

#Filter out classes where there are less than 10 objects (for density plot)
small_body_dt <- small_body_dt[, .SD[.N>10], by=class]

ggplot(small_body_dt, aes(x=moid)) + 
  geom_density(aes(color=class)) +
  labs(x='Minimum distance from Earth (au)', y='Density distribution (arb)')

ggplot(small_body_dt, aes(x=diameter)) + 
  geom_density(aes(color=class)) + coord_cartesian(xlim = c(0, 10)) +
  labs(x='Diameter (km)', y='Density distribution (arb)')
  
############################################################################
#Online follow up
  
#Kepler's 3rd log-log plot
ggplot(test, aes(x = x**2, y = y**3)) +
  geom_point() +
  xlab("a [AU]") + ylab("Period [year]")

ggplot(small_body_dt[51:100], aes(x = per_y**2, y = a**3)) +
            geom_point() +
            scale_x_log10() + scale_y_log10() +
            xlab("a [AU]") + ylab("Period [year]")
  
#Fit of kepler's 3rd
lm(log10(per_y) ~ log10(a), data =small_body_dt)

ggplot(small_body_dt, aes(x = e)) +
     geom_histogram(binwidth = 0.05, fill = "lightblue", colour = "black") +
     xlab("Eccentricity") + ylab("Count") +
     theme_classic()

ggplot(small_body_dt, aes(x = e, y = a)) +
     geom_point(aes(color = class), alpha = 0.5, size = 2) +
     scale_y_log10() +
     scale_colour_discrete("Class", h = c(240, 0), c = 150) +
     xlab("Eccentricity") + ylab("Semi-Major Axis [AU]") +
     theme_classic()

library(nleqslv)
library(plyr)

small_body_dt <- ddply(small_body_dt, .(full_name), mutate, E = {
     kepler <- function(E) {
         ma - E + e * sin(E)
       }
       nleqslv(0, kepler, method = "Newton")$x
   })

#install.packages('scales')
library(scales)

degrees <- function(x, ...) {parse(text = paste0(x, "*degree"))}

ggplot(small_body_dt[101:200,], aes(x = ma / pi * 180, y = E / pi * 180, colour = e)) +
  geom_point() +
  scale_colour_gradient("Eccentricity (e)", limits=c(0, 1), high = "red", low = "blue") +
  xlab("Mean Anomaly (M)") + ylab("Eccentric Anomaly (E)") +
  theme_classic()

small_body_dt <- transform(small_body_dt, theta = 2 * atan2(sqrt(1 + e) * sin (E/2), sqrt(1 - e) * cos (E/2)))

small_body_dt <- within(small_body_dt, {Z = 0 
  Y = a * sqrt(1 - e**2) * sin(E) 
  X = a * cos(E) - a * e})

transformation.matrix <- function(w, om, i) {
  a11 =   cos(w) * cos(om) - sin(w) * cos(i) * sin(om)
  a12 =   cos(w) * sin(om) + sin(w) * cos(i) * cos(om)
  a13 =   sin(w) * sin(i)
  a21 = - sin(w) * cos(om) - cos(w) * cos(i) * sin(om)
  a22 = - sin(w) * sin(om) + cos(w) * cos(i) * cos(om)
  a23 =   cos(w) * sin(i)
  a31 =   sin(i) * sin(om)
  a32 = - sin(i) * cos(om)
  a33 =   cos(i)
  matrix(c(a11, a12, a13, a21, a22, a23, a31, a32, a33), nrow = 3, byrow = TRUE)
  }

#Try for one object
test <- (A = do.call(transformation.matrix, small_body_dt[1, c("w", "om", "i")]))

#Apply to all objects
small_body_dt <-  ddply(small_body_dt, .(full_name), function(df) {
A = transformation.matrix(df$w, df$om, df$i)

# Invert transformation matrix
A = t(A)
# Transform to reference frame
r = A %*% matrix(unlist(df[, c("X", "Y", "Z")]), ncol = 1)
r = matrix(r, nrow = 1)
colnames(r) = c("x", "y", "z")
cbind(df, r)
})

small_body_dt %<>%
  mutate(radius = (x**2 + y**2 + z**2)**0.5)

library(grid)
#Dashed line shows 1 AU radius
#Sampling only moid < 0.001
ggplot(small_body_dt[small_body_dt$moid<0.1 & small_body_dt$moid>0.09,], aes(x = x, y = y)) +
  geom_point(alpha=0.25) +
  xlab("x [AU]") + ylab("y [AU]") +
  scale_x_continuous(limits = c(-5, 5)) + scale_y_continuous(limits = c(-5, 5)) +
  annotation_custom(grob = circleGrob(
  r = unit(0.5, "npc"), gp = gpar(fill = "transparent", col = "black", lwd = 2, lty = "dashed")),
  xmin=-1, xmax=1, ymin=-1, ymax=1) +
  theme_classic() +
  theme(aspect.ratio = 1)

ggplot(small_body_dt[small_body_dt$moid<0.001,], aes(x = x, y = z)) +
  geom_point(alpha = 0.25) +
  xlab("x [AU]") + ylab("z [AU]") +
  scale_x_continuous(limits = c(-5, 5)) + scale_y_continuous(limits = c(-5, 5)) +
  theme_classic() +
  theme(aspect.ratio = 1)

#3d chart
#install.packages("rgl")
library(rgl)

temp <- small_body_dt[small_body_dt$moid<0.1,]
class_temp <- unique(temp$class)
col_temp <- heat.colors(length(class_temp), alpha=NULL)
class_col = c(class_temp=col_temp)

plot3d(temp[abs(temp$radius) < 2,'x'],
       temp[abs(temp$radius) < 2,'y'],
       temp[abs(temp$radius) < 2,'z'],
       col = class_col) + 
  spheres3d(x=0, radius = 1, alpha=0.5, color='grey') +
  spheres3d(x=0, radius = 0.2, alpha = 1, color='green')

class_temp <- unique(small_body_dt$class)
col_temp <- heat.colors(length(class_temp), alpha=NULL)
class_col = c(class_temp=col_temp)

plot3d(small_body_join[,'x'],
       small_body_join[,'y'],
       small_body_join[,'z'],
       col = class_col) + 
  spheres3d(x=0, radius = 1, alpha=0.5, color='blue') +
  spheres3d(x=0, radius = 0.1, alpha = 1, color='green')

plot3d(x=unlist(small_body_dt[small_body_dt$radius < 10,'x']),
       y=unlist(small_body_dt[small_body_dt$radius < 10,'y']),
       z=unlist(small_body_dt[small_body_dt$radius < 10,'z'])) 
  #spheres3d(x=0, radius = 1, alpha=0.5, color='blue') +
  #spheres3d(x=0, radius = 0.1, alpha = 1, color='green')

#Schematic
# 1. 3D plot (one or two times, static)
# 2. Interactive histograms/barplots for EDA
# 3. Damage estimator with probability

small_body_join <- small_body_join[, .SD[.N>10], by=class]

small_body_join <- small_body_join %>% 
  mutate(., round_dradius = round(radius*2))

small_body_join <- small_body_join %>% 
  filter(radius < 50)

small_body_dt <- small_body_dt %>% 
  filter(radius < 50)

###################################################################
#First try at heat map

#These need to be reactive based on radius filtering
minmax_x <- small_body_join %>% 
  select(x) %>% 
  summarise(min(.),
            max(.))
seq_x <- seq(minmax_x[[1]]/2, minmax_x[[2]]/2, length.out = 11)
#This also needs to be based on radius filtering
minmax_y <- small_body_join %>% 
  select(y) %>% 
  summarise(min(.),
            max(.))
seq_y <- seq(minmax_y[[1]], minmax_y[[2]], length.out = 11)

in_interval <- function(x, interval){ 
  stopifnot(length(interval) == 2L) 
  interval[1] < x & x < interval[2] 
} 

find_location <- function(x, seq_x){
  ifelse(in_interval(x, c(seq_x[1], seq_x[2])), 1,
         ifelse(in_interval(x, c(seq_x[2], seq_x[3])), 2,
         ifelse(in_interval(x, c(seq_x[3], seq_x[4])), 3,
         ifelse(in_interval(x, c(seq_x[4], seq_x[5])), 4,
         ifelse(in_interval(x, c(seq_x[5], seq_x[6])), 5,
         ifelse(in_interval(x, c(seq_x[6], seq_x[7])), 6,
         ifelse(in_interval(x, c(seq_x[7], seq_x[8])), 7,
         ifelse(in_interval(x, c(seq_x[8], seq_x[9])), 8,
         ifelse(in_interval(x, c(seq_x[9], seq_x[10])), 9,
         10)))))))))
}
test <- small_body_dt %>% 
  filter(class == 'MBA')
x_loc <- find_location(test$x, seq_x)
y_loc <- find_location(test$y, seq_y)


mat_loc <- matrix(0L, nrow = 10, ncol = 10)
mat_loc <- data.frame(x_loc, y_loc)

mat_loc <- mat_loc %>% 
  group_by(x_loc, y_loc) %>% 
  dplyr::summarise(count = n())

ggplot(data = mat_loc, aes(x = x_loc, y = y_loc, fill = count)) + geom_tile()

#########################################################################
#Second try at heatmap
install.packages('akima')
require(akima)

data <- data.frame(x=small_body_join$x,
                   y=small_body_join$y,
                   distance=small_body_join$moid)

resolution <- 0.1

a <- interp(x=data$x, y=data$y, z=data$distance, 
            xo=seq(min(data$x),max(data$x),by=resolution), 
            yo=seq(min(data$y),max(data$y),by=resolution), duplicate="mean")

filled.contour(a, color.palette=heat.colors)






#Select class of meteor
#Display heatmap/plot3D cached image of orbit
#Blurb for class description
#
#Compare/contrast all asteroids via e, a, diameter...

####MBA = middle belt asteroid
#OMB = outer belt asteroid
#MCA = mars crossing asteroid
#APO = Apollo NEOs
#AMO = Amor NEOs
#IMB = inner belt asteroid
#ATE = Aten NEOs
#TNO = trans-neptunian objects
#CEN = centaur
#TJN = jupiter trojan
#IEO = Atira NEOs
#AST = "other" asteroid
#COM = "other" comet
#HTC = Halley-type comets
#JFc = Jupiter family comets (old)
#JFC = Jupiter family comets (new)
#CTc = chiron type
#ETc = Enke type comets

meteor_classes <- unique(small_body_dt$class)
meteor_blurb <- c(paste(tags$b("Outer belt asteroid"),
                        tags$br(),
                        'Lie within Zone III of the Kirkwood schema.
                        Semi-major axes lie beyond 2.82 AU.'),
                  paste(tags$b("Mars crossing asteroid: "),
                        tags$br(),
                        'Asteroids whose orbit crosses that of Mars. 
                        Many databases, only list asteroids with a perihelion greater than
                        1.3 AU as Mars-crossers. An asteroid with a perihelion less than this
                        is classed as a near-Earth object even though it is crossing the orbit
                        of Mars as well as crossing (or coming near to) that of Earth.'),
                  paste(tags$b("Apollo asteroids: "),
                        tags$br(),
                        'A group of near-Earth asteroids named after 1862 Apollo,
                        discovered by German astronomer Karl Reinmuth in the 1930s. They are Earth crossing 
                        asteroids that have an orbital semi-major axis greater than that of the Earth (> 1 AU)
                        but perihelion distances less than the Earth\'s aphelion distance (q < 1.017 AU).'),
                  paste(tags$b("Amor type object: "),
                        tags$br(),
                        'A group of near-Earth asteroids named after the asteroid
                        1221 Amor. They approach the orbit of Earth from beyond, but do not cross it. 
                        Most Amors cross the orbit of Mars.',
                        tags$br(),
                        'a > 1.0 AU and 1.017 AU < q < 1.3 AU'),
                  paste(tags$b("Inner belt asteroid: "),
                        tags$br(),
                        'Lie within Zone I of the Kirkwood schema.
                        Semi-major axes lie between from 2.06 to 2.5 AU'),
                  paste(tags$b("Aten type object: "),
                        tags$br(),
                        'A group of near-Earth asteroid with a semi-major axis of less than one astronomical unit (AU). 
                        They have an aphelion (furthest distance from the Sun) greater than 0.983 AU. 
                        An Aten orbit need not be entirely contained within Earth\'s orbit; nearly all known Aten 
                        asteroids have aphelion greater than 1 AU and semi-major axis less than 1 AU.'),
                  paste(tags$b("Trans-Neptunian objects: "),
                        tags$br(),
                        'Any minor planet in the Solar System that orbits the Sun at a greater average distance (semi-major axis)
                        than Neptune, 30 astronomical units (AU). Twelve minor planets with a semi-major axis greater than 150 AU
                        and perihelion greater than 30 AU are known, which are called extreme trans-Neptunian objects (ETNOs).'),
                  paste(tags$b("Centaur type object: "),
                        tags$br(),
                        'Small bodies with size similar to asteroids but 
                        composition similar to comets, that revolve around
                        the Sun in the outer solar system, 
                        mainly between the orbits of Jupiter and Neptune.'),
                  paste(tags$b("Jupiter Trojan object: "),
                        tags$br(),
                        'Large group of asteroids that share the orbit of the planet Jupiter around the Sun.
                        Relative to Jupiter, each Trojan librates around one of Jupiter\'s two stable Lagrangian
                        points, L4, lying 60 degrees ahead of the planet in its orbit, and L5, 60 degrees behind.
                        Jupiter trojans are distributed in two elongated, curved regions around these Lagrangian
                        points with an average semi-major axis of about 5.2 AU.'),
                  paste(tags$b("Atira type object: "),
                        tags$br(),
                        'Also known as Interior-Earth Objects (IEOs), asteroids whose orbits are entirely
                        confined within Earth\'s orbit, that is, their orbit has an aphelion (farthest point
                        from the Sun) smaller than Earth\'s perihelion (nearest point to the Sun),
                        which is 0.983 astronomical units (AU). '),
                  paste(tags$b("Undefined asteroid: "),
                        tags$br(),
                        'Characteristics unknown, yet to be determined due to lack of observations.'),
                  paste(tags$b("Undefined comet: "),
                        tags$br(),
                        'Characteristics unknown, yet to be determined due to lack of observations.'),
                  paste(tags$b("Halley-type comet: "),
                        tags$br(),
                        'Have orbital periods between 20 and 200 years, and orbits which can be highly
                        inclined to the ecliptic. Both of these observations suggest that Halley-type
                        comets have a different origin to the other class of short-period comet, the
                        Jupiter-family comets, which have shorter periods and lower inclinations.'),
                  paste(tags$b("Jupiter family comet (old): "),
                        tags$br(),
                        'Named as such since their current orbits are primarily determined by the
                        gravitational influence of Jupiter. Although the orbits of Jupiter-family
                        comets are contained within or do not extend much beyond the orbit of Jupiter,
                        it is thought that they originate in the Kuiper Belt, a collection of rock-ice
                        bodies located just beyond the orbit of Neptune. Orbital period < 20y'),
                  paste(tags$b("Jupiter family comet (new): "),
                        tags$br(),
                        'Named as such since their current orbits are primarily determined by the
                        gravitational influence of Jupiter. Although the orbits of Jupiter-family
                        comets are contained within or do not extend much beyond the orbit of Jupiter,
                        it is thought that they originate in the Kuiper Belt, a collection of rock-ice
                        bodies located just beyond the orbit of Neptune. 2 < T_jup < 3'),
                  paste(tags$b("Chiron type comet: "),
                        tags$br(),
                        'T_j > 3 and a > a_j'),
                  paste(tags$b("Enke type comet: "),
                        tags$br(),
                        'T_j > 3 and a < a_j'))


meteor_descriptions <- data.frame(meteor_classes, meteor_blurb, stringsAsFactors = F)
fwrite(meteor_descriptions, 'meteor_descriptions')

mba_dt <- small_body_dt %>% filter(class == 'MBA')
small_body_dt <- small_body_dt %>% filter(class != 'MBA')





###############################################################
#Crater equation
#From: 'http://keith.aa.washington.edu/craterdata/scaling/theory.pdf'

#V_cr = crater volume
#rho_t = target mass density (g/cm^3)
#mass = mass
## pi_v = rho_t*V_cr/mass

#k_1 = coeff_1
#k_2 = coeff_2
#y_t = target strength (dynes/cm^3)
#u_s = source velocity (*NEEDS to be in cm/s*)
#u_f = u_s*cos(theta) = final velocity
#a_s = source radius (*NEEDS to be in cm*)
#delta_s = source mass density (g/cm^3)
#mu and nu are coefficients of equation


crater_formation <- function(a_s, u_s, rho_t, delta_s, y_t, mu, nu, k_1, k_2, k_r, k_d){
  g = 980.7 #cm/s^2
  #Coefficient calc pi_2
  pi_2 <- g*a_s/u_s**2

  #Coefficient calc pi_3
  pi_3 <- y_t/(rho_t*u_s**2)

  #Coefficient calc pi_v
  pi_v <- k_1*(pi_2(g, a_s, u_s)*(rho_t/delta_s)**((6*nu-2-mu)/(3*mu)) +
               (k_2*pi_3(y_t, rho_t, u_s)*(rho_t/delta_s)**((6*nu-2)/(3*mu)))**((2+mu)/2))**(-3*mu/(2+mu))

  #Crater volume V_cr
  V_cr <- pi_v*(pi*(4/3)*(a_s)**3 * delta_s)/rho_t

  #Crater radius r_cr
  r_cr <- k_r*V_cr**(1/3)

  #Crater depth d_cr
  d_cr <- k_d*V_cr**(1/3)

  V_ej = 0.8*V_cr #ejected volume
  
  T_form = 0.8*(V_cr**(1/3)/g)**0.5 #Formation time
  
  value_cm <- c(V_cr, V_ej, r_cr, d_cr, T_form)
  value_km <- value_cm*c((1e-5)**3, (1e-5)**3, 1e-5, 1e-5, 1)
  value_mi <- value_km*c(0.6124**3, 0.6124**3, 0.6124, 0.6124, 1)
  
  output <- data.frame(name=c('V_cr', 'V_ej', 'r_cr', 'd_cr', 'T_form'),
                       value_cm=value_cm,
                       value_km=value_km,
                       value_mi=value_mi,
                       stringsAsFactors = F)
  
  return(output)
}

materials_name <- c('water',
                    'dry sand',
                    'soil',
                    'soft rock',
                    'hard rock',
                    'cold ice')

materials_k_1 <- c(0.904,
                   0.14,
                   0.094,
                   0.012,
                   0.012,
                   0.012)

materials_k_2 <- c(0,
                   0,
                   0.75,
                   0.8,
                   0.8,
                   0.8)

materials_mu <- c(0.55,
                   0.4,
                   0.55,
                   0.55,
                   0.55,
                   0.55)

materials_nu <- c(0.33,
                   0.33,
                   0.33,
                   0.33,
                   0.33,
                   0.33)

materials_y_t <- c(0,
                   0,
                   1.2e6,
                   1.5e8,
                   1.6e9,
                   1.5e8)

materials_rho_t <- c(1,
                   1.7,
                   2.1,
                   2.1,
                   3.2,
                   0.93)

materials_k_r <- c(0.8,
                  1.4,
                  1.1,
                  1.1,
                  1.1,
                  1.1)

materials_k_d <- c(0.75,
                   0.35,
                   0.6,
                   0.6,
                   0.6,
                   0.6)

materials <- data.frame(name = materials_name,
                        k_1 = materials_k_1,
                        k_2 = materials_k_2,
                        mu = materials_mu,
                        nu = materials_nu,
                        rho_t = materials_rho_t,
                        y_t = materials_y_t,
                        k_r = material_k_r,
                        k_d = materials_k_d)


impactor_name <- c('aluminum (2.7)',
                   'plastic (0.95)',
                   'steel (7.8)',
                   'c-type (1.8)',
                   's-type (3.0)',
                   'comet (0.8)')

impactor_delta_s <- c(2.7,
                      0.95,
                      7.8,
                      1.8,
                      3.0,
                      0.8)

impactor <- data.frame(name = impactor_name,
                       delta_s = impactor_delta_s)




target <- materials %>% filter(materials_name == 'water')
crater <- crater_formation(k_1 = target$k_1,
                         k_2 = target$k_2,
                         mu = target$mu,
                         nu = target$nu,
                         rho_t = target$rho_t,
                         y_t = target$y_t,
                         k_r = target$k_r,
                         k_d = target$k_d,
                         u_s = 1e5*meteor_vel, #Meteor velocity in km/s trans to cm/s
                         a_s = 0.5e5*meteor_diam, #Meteor diameter in km trans to radius of cm
                         delta_s = impactor$delta_s[meteor_mat])

m <-leaflet() %>%
  addTiles() %>% 
  fitBounds(lng-0.5, lat-0.5, lng+0.5, lat+0.5) %>% 
  addCircles(lng = lng, lat = lat, radius = meteor_diam) %>% 
  addPopups(lng = lng, lat = lat, popup=paste(city_name, '-- Crater diameter (km) = ', meteor_diam))
  
city_name <- c("New York",
               "London",
               "Moscow",
               "Beijing",
               "Sydney")

city_lat <- c(40.712784,
              51.507351,
              55.755826,
              39.904211,
              -33.868820)

city_lon <- c(-74.005941,
              -0.127758,
              37.617300,
              116.407395,
              151.209296)

city_dt <- data.frame(name=city_name,
                      lat=city_lat,
                      lon=city_lon)

fwrite(sentry_dt, 'sentry_dt.dt')
fwrite(small_body_dt, 'small_body_dt.dt')
fwrite(small_body_join, 'small_body_join.dt')
fwrite(sbdt_summary, 'sbdt_summary.dt')
fwrite(mba_dt, 'mba_dt.dt')

fwrite(impactor, 'impactor')
fwrite(materials, 'materials')
fwrite(city_dt, 'city_dt.dt')


ggplot(data=mba_dt, aes(x=a)) +
  geom_density(fill='yellow', alpha=0.2)


df_sun <- ellipse_create(a = 0.1, q = 0.1) %>% mutate(out2.html.tooltip = 'Sun')
df_earth <- ellipse_create(a = 1, q = 1) %>% mutate(out2.html.tooltip = 'Earth')
df_mars <- ellipse_create(a = 1.3813, q = 1.5237) %>% mutate(planet = 'Mars')
df_jupiter <- ellipse_create(a = 4.952, q = 5.203) %>% mutate(planet = 'Jupiter')

sun_plot <- gvisScatterChart(df_sun,
                             options=list(lineWidth=2, pointSize=0, colors="['yellow']"))
earth_plot <- gvisScatterChart(df_earth,
                               options=list(lineWidth=2, pointSize=0, colors="['blue']"))
plot(gvisMerge(sun_plot,earth_plot))


AMO <- sbdt_summary[1]
APO <- sbdt_summary[2]
AST<- sbdt_summary[3]
ATE<- sbdt_summary[4]
CEN<- sbdt_summary[5]
COM<- sbdt_summary[6]
CTc<- sbdt_summary[7]
ETc<- sbdt_summary[8]
HTC<- sbdt_summary[9]
IEO<- sbdt_summary[10]
IMB<- sbdt_summary[11]
JFc<- sbdt_summary[12]
JFC<- sbdt_summary[13]
MBA<- sbdt_summary[14]
MCA<- sbdt_summary[15]
OMB<- sbdt_summary[16]
TJN<- sbdt_summary[17]
TNO<- sbdt_summary[18]


orbital_plot + 
  geom_path(data = ellipse_create(a = AMO$avg_a,
                                  q = AMO$avg_q),
            color = 'black') +
  geom_path(data = ellipse_create(a = APO$avg_a,
                                  q = APO$avg_q),
            color = 'black') +
  #geom_path(data = ellipse_create(a = AST$avg_a,
  #                                q = AST$avg_q),
  #          color = 'black') +
  geom_path(data = ellipse_create(a = ATE$avg_a,
                                  q = ATE$avg_q),
            color = 'black') +
  #geom_path(data = ellipse_create(a = CEN$avg_a,
  #                                q = CEN$avg_q),
  #          color = 'black') +
  #geom_path(data = ellipse_create(a = COM$avg_a,
  #                                q = COM$avg_q),
  #          color = 'black') +
  #geom_path(data = ellipse_create(a = CTc$avg_a,
  #                                q = CTc$avg_q),
  #          color = 'black') +
  #geom_path(data = ellipse_create(a = ETc$avg_a,
  #                                q = ETc$avg_q),
  #          color = 'black') +
  #geom_path(data = ellipse_create(a = HTC$avg_a,
  #                                q = HTC$avg_q),
  #          color = 'black') +
  geom_path(data = ellipse_create(a = IEO$avg_a,
                                  q = IEO$avg_q),
            color = 'black') #+
  #geom_path(data = ellipse_create(a = IMB$avg_a,
  #                                q = IMB$avg_q),
  #          color = 'black') +
  #geom_path(data = ellipse_create(a = JFc$avg_a,
  #                                q = JFc$avg_q),
  #          color = 'black') +
  #geom_path(data = ellipse_create(a = JFC$avg_a,
  #                                q = JFC$avg_q),
  #          color = 'black') +
  #geom_path(data = ellipse_create(a = MBA$avg_a,
  #                                q = MBA$avg_q),
  #          color = 'black') +
  #geom_path(data = ellipse_create(a = MCA$avg_a,
  #                                q = MCA$avg_q),
  #          color = 'black') +
  #geom_path(data = ellipse_create(a = OMB$avg_a,
  #                                q = OMB$avg_q),
  #          color = 'black') +
  #geom_path(data = ellipse_create(a = TJN$avg_a,
  #                                q = TJN$avg_q),
  #          color = 'black') #+
  #geom_path(data = ellipse_create(a = TNO$avg_a,
  #                                q = TNO$avg_q),
  #          color = 'black')
  
ggplotly(orbital_plot + 
           geom_path(data = ellipse_create(a = AMO$avg_a,
                                           q = AMO$avg_q),
                     color = 'black')) %>% 
  add_annotations(x = ellipse_create(a = AMO$avg_a,
                                     q = AMO$avg_q)$out1[75],
                  y = ellipse_create(a = AMO$avg_a,
                                     q = AMO$avg_q)$out2[75],
                  text = 'test',
                  xref = "x",
                  yref = "y",
                  showarrow = TRUE,
                  arrowhead = 4,
                  arrowsize = .5,
                  ax = 20,
                  ay = -40)


add_annotations(x = ellipse_create(a = temp$avg_a,
                                   q = temp$avg_q)$out1[75],
                y = ellipse_create(a = temp$avg_a,
                                   q = temp$avg_q)$out2[75],
                text = paste(temp$class, ' orbit'),
                xref = "x",
                yref = "y",
                showarrow = F,
                bgcolor = 'white',
                bordercolor = 'black')
