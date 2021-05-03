# This is the first script

df <- load("some_data.Rdata")

fit <- lm(some ~ new_analyses,data=df)

summary(fit)

plot(fit)
