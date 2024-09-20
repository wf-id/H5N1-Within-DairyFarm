library(deSolve)
library(ggplot2)
library(dplyr)

##################################
#Set initial conditions and params
##################################
init.inf = 1 #initial infected
herd.size = 100 #total herd size
r0 = 1.2 #R0
fever.duration = 1.8 #infectious period
milk.loss.dur = 7 #recovery 1 period
cow.prod.lifespan = 365.25*2 #lifespan for popn turnover
include.prod.lifespan = FALSE #initial can toggle on or off in shiny
base.milk.prod = 100 # in units of lbs per day
sympt.milk.prod = 100-11.4 #reduced milk production
discard.sick.prod = TRUE #can toggle on and off in shiny
disease.induced.mortality = .1 #death rate due to H5N1
prop.cull = 0.1 #death rate due to culling
sim.length = 365 #simulation length
add_milk = FALSE
state <- c("S" = 99, # Symptomatic compartment
           "I" = 1, # Infectious compartment
           "B" = 0, # Symptomatic, but not infectious compartment
           "R" = 0, # Recovered compartment
           "D" = 0, # Death from disease
           "C" = 0, # Culled compartment
           "Z" = 0, # Dummy compartment to count cows that were infected but died due to productive lifespan
           "M" = 0)  # Milk compartment

# Set time in units of days
init.time <- 0
end.time <- 365
t.step.dur <- 1 / 24 # time step of an hour
time <- seq(from = init.time, to = end.time, by = t.step.dur)

# Cow parameters
parms <- c(beta = r0 / (fever.duration) / herd.size,
           gamma = 1 / fever.duration,
           alpha = 1 / (milk.loss.dur - fever.duration),
           mu = (1 / cow.prod.lifespan) * dplyr::if_else(include.prod.lifespan, 1, 0),
           phi = disease.induced.mortality,
           m = base.milk.prod,
           p = 1 - sympt.milk.prod / base.milk.prod,
           c = prop.cull)

##################################
#SIR model in absolute time
##################################
SIRMmod <- function (Time, State, pars) {
  with(as.list(c(State, pars)), {
    N <- sum(State) - M - C - D - Z
    #I <- Ip + Is
    dS = -beta * S * I + N * mu - S * mu
    #! TODO remove this compartment
    dI = beta * S * I - gamma * I - I * mu
    dB = (gamma * I) * (1 - phi) - alpha * B - mu * B
    dR =  (1 - c) * alpha * B - mu * R
    dC = alpha * c * B
    dD = gamma * I * phi
    dZ = mu*R + mu*I + mu*B
    dM = m * ((S + R + I) + (B) * (1 - p))
    return(list(c(dS, dI, dB, dR, dD, dC, dZ, dM)))
  })
}

#run model
out <- ode(state, time, SIRMmod, parms, rtol = 1e-15, maxsteps = 500000)
sub <- as.data.frame(out)%>%
  filter(time < 30)

##################################
#Figure
##################################
Pal1 <- c('Susceptible' = 'dodgerblue3',
          'Infectious' = 'red3',
          'Clinical, not infectious' = 'orange1',
          'Recovered' = 'gold',
          'Died' = 'black',
          'To market'= 'gray50')

fig <- ggplot()+
  theme_classic(base_size = 36)+ #white background w/ out grid
  geom_line(data=sub, aes(x=time, y=S, color="Susceptible"), size=1.5, alpha=0.8)+
  geom_line(data=sub, aes(x=time, y=I, color="Infectious"), size=1.5, alpha=0.8)+
  geom_line(data=sub, aes(x=time, y=B, color="Clinical, not infectious"), size=1.5, alpha=0.8)+
  geom_line(data=sub, aes(x=time, y=R, color="Recovered"), size=1.5, alpha=0.8)+
  geom_line(data=sub, aes(x=time, y=D, color="Died"), size=1.5, alpha=0.8)+
  geom_line(data=sub, aes(x=time, y=C, color="To market"), size =1.5, alpha=0.8)+
  scale_color_manual(values = Pal1, name= "Disease state",
                     breaks=c("Susceptible","Infectious","Clinical, not infectious", "Recovered", 'Died', 'To market')) +
  # theme(text = element_text(size=20),
  #       axis.text.x = element_text(size=16))+
  labs(y= "% of herd", x = "Time (days)") #axis labels

ggsave("ModelOutput.tiff", plot=fig, path = here::here("figures"), width = 10, height = 4, device='tiff', dpi=300)

