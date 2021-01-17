rm(list = ls())

require("ggplot2"); require("gridExtra"); require("ks")

cv <- function(x){sd(x)/mean(x)*100}

data("potthoffroy", package = "mice")
dados <- potthoffroy; rm(potthoffroy)

dados <- data.frame(id = rep(dados$id, times = 4), sex = rep(dados$sex, times = 4),
                    age = rep(c(8, 10, 12, 14), each = 27), y = vec(dados[,3:6]))
dados$medias = rep(tapply(dados$y, dados$age, mean), table(dados$age))
dados$dps = rep(tapply(dados$y, dados$age, sd), table(dados$age))

dadosF <- dados[dados$sex == "F",]
dadosF$medias = rep(tapply(dadosF$y, dadosF$age, mean), table(dadosF$age))
dadosF$dps = rep(tapply(dadosF$y, dadosF$age, sd), table(dadosF$age))

dadosM <- dados[dados$sex == "M",]
dadosM$medias = rep(tapply(dadosM$y, dadosM$age, mean), table(dadosM$age))
dadosM$dps = rep(tapply(dadosM$y, dadosM$age, sd), table(dadosM$age))


plot0 <- ggplot(dados, aes(x = age, y = y, group = id, colour = sex)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ylab("Distância glândula pituitária da fissura pterigomaxilar (mm)") +
  xlab("Idade (anos)") +
  scale_color_manual(values = c("indianred1", "lightseagreen"), labels =c("Feminino","Masculino")) +
  labs(color='Sexo') 
plot0
  

plot1 <- ggplot(dados, aes(x = age, y = y, group = id, colour = sex)) +
  geom_line() +
  geom_point() +
  theme_bw() +
  ylab("Distância glândula pituitária da fissura pterigomaxilar (mm)") +
  xlab("Idade (anos)") +
  scale_color_manual(values = c("indianred1", "lightseagreen", "black"), labels =c("Feminino","Masculino", "Perfil Medio")) +
  labs(color='Sexo') +
  geom_line(aes(y = medias, colour = "Perfil Medio"), size = 0.8)

p10 <- plot0 + facet_wrap(~sex, ncol = 2)
p11 <- plot1 + facet_wrap(~sex, ncol = 2)

box <- ggplot(dados, aes(x = factor(age), y = y, fill = sex)) +
  geom_boxplot() +
  scale_fill_manual(values = c("indianred1", "lightseagreen"), labels =c("Feminino","Masculino")) +
  labs(fill='Sexo') +
  xlab("Idade (anos)") +
  theme_bw()

p2 <- ggplot(dados, aes(x=age, y=medias)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=medias-dps, ymax=medias+dps), width=.2,
                position=position_dodge(0.05)) +
  ylab("Distância glândula pituitária da fissura pterigomaxilar (mm)") +
  xlab("Idade (anos)") +
  theme_bw(); p2

p21 <- ggplot(dados, aes(x=age, y=medias)) + 
  geom_line() +
  geom_point()+
  geom_errorbar(aes(ymin=medias-dps, ymax=medias+dps), width=.2,
                position=position_dodge(0.05)) +
  ylab("Distância glândula pituitária da fissura pterigomaxilar (mm)") +
  xlab("Idade (anos)") +
  facet_wrap(~sex, ncol = 2) +
  scale_fill_manual(values = c("indianred1", "lightseagreen")) +
  scale_color_manual(values = c("indianred1", "lightseagreen")) +
  theme_bw(); p21

