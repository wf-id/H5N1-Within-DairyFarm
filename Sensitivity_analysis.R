library(dplyr)
library(ggplot2)

ffi.obs=seq(0.3,0.8, 0.1) #sequence
N = 500

r0_se <- function(r, z, N){
  se <- sqrt((1 + (1-z)* r^2)/ (N * z * (1-z)))
  data.frame(estimate = round(r,2), z, ci_lo  = round(r - 2 * se, 2), ci_high = round(r + 2* se,2))
  #gives estimate of standard error from stochastic models
  #r = sample mean of R0, z = FFI, N =  herd size N - confirm not # of stochastic iterations
}

s_inf <- function(ffi_e, prop_susceptible){
  (1-ffi_e)*prop_susceptible
}

df <- data.frame(prop_susceptible = rep(seq(min(ffi.obs), 1, 0.01), length(ffi.obs)))%>% #proportion of the herd susceptible
#  mutate(ffi_obs = c(rep(ffi.low, (length(prop_susceptible))/2),rep(ffi.high, (length(prop_susceptible))/2)), #observed FFI
   mutate(ffi_obs = c(rep(ffi.obs, each=length(prop_susceptible)/length(ffi.obs))), #observed FFI
         ffi_e = ffi_obs/prop_susceptible, #effective FFI
         ffi_e = round(ffi_e, 2),
         R0_estimate = -log(s_inf(ffi_e, prop_susceptible))/(1-s_inf(ffi_e, prop_susceptible)),
         R0_estimate = round(R0_estimate, 2))

df.se <- r0_se(df$R0_estimate, df$ffi_e, N*df$prop_susceptible) %>% rename("R0_estimate" = "estimate", "ffi_e" = "z")
df<- full_join(df, df.se)

#sensitivity analysis full range
fig <- ggplot()+
  theme_classic(base_size = 24)+
  geom_line(data=df,  aes(x= prop_susceptible, y= R0_estimate, color=as.character(ffi_obs)), size=2)+
  scale_y_continuous(limits=c(0.9,5))+
  geom_ribbon(data=df, aes(x=prop_susceptible, ymin=ci_lo,ymax=ci_high, fill= as.character(ffi_obs)), alpha=0.5)+
  scale_colour_viridis_d("Observed\nFFI")+ scale_fill_viridis_d("Observed\nFFI")+
  xlab("Proportion susceptible") + ylab(bquote(''*R[0]*' estimate'))


#sensitivity analysis just 0.4 FFI
fig1 <- ggplot()+
  theme_classic(base_size = 24)+
  geom_vline(xintercept = 0.4, col = "black", linetype=3, size=2)+
  geom_line(data=df %>% filter(ffi_obs==0.4),  aes(x= prop_susceptible, y= R0_estimate, color=as.character(ffi_obs)), size=2)+
  #scale_y_continuous(limits=c(0.9,5))+
  scale_y_continuous(breaks = c(0, 1.2, 3, 5), limits = c(0, 5)) +

  geom_ribbon(data=df %>% filter(ffi_obs==0.4), aes(x=prop_susceptible, ymin=ci_lo,ymax=ci_high, fill= as.character(ffi_obs)), alpha=0.5)+
  geom_hline(yintercept = 1.2, col = "red3", linetype=2, size=2)+
  scale_colour_viridis_d("Observed\nFFI")+ scale_fill_viridis_d("Observed\nFFI")+
  xlab("Proportion susceptible") + ylab(bquote(''*R[0]*' estimate'))


#RO among different groups
df <- data.frame(Group= c("None", "Dry", "1st", "2nd", "3+"),
                 R0 = c(0.25, 0.75, 2.5, 3.5, 4))
fig2 <- ggplot()+
  theme_classic(base_size=24)+
  scale_y_continuous(breaks = c(0,1.2, 5), limits = c(0, 5), labels= c("low", "1.2", "high")) +
  geom_hline(yintercept = 1.2, col = "red3", linetype=2, size=2)+
  geom_point(data=df, aes(y=R0, x=factor(Group, Group)), fill=NA, shape=22, size=20)+
  ylab(bquote(''*R[0]*' estimate'))+
  xlab('Lactation')


fig <- ggpubr::ggarrange(fig2, fig1, ncol=2,
                  widths=c(8, 10),
                  labels = c("A", "B"))

ggsave("R0_sensitivity.tiff", plot=fig, path = here::here("figures"), width = 10, height = 4, device='tiff', dpi=300)
