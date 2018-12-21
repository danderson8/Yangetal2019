library(lmtest)
library(qvalue)

data <- read.csv("PTM_converted.csv")
attach(data)

sink(file="ROutput.txt")

write("Terms, PTM_Zn", file="PTM_Terms.csv")
write("TermsR2, PTM_ZnR2", file="PTM_R2_Terms.csv")

act[which(act <= 0.0001)] = 0.0001

DHC <- log(act, base=10) # /(0.58*1300)*1000000/60

lm_pos1 <- lm(DHC ~ pos72)

print("Position 72")
print(coefficients(lm_pos1))
print(summary(lm_pos1)$coefficients)
print(summary(lm_pos1))
write(c("pos72",coefficients(lm_pos1)[2]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
print(summary(lm_pos1)$r.squared)
write(c("pos72",summary(lm_pos1)$r.squared), sep=",", file="PTM_R2_Terms.csv", append=TRUE, ncol=2)

p_value_pos1 <- summary(lm_pos1)$coefficients[2,4]


lm_pos2 <- lm(DHC ~ pos193)

print("Position 193")
print(coefficients(lm_pos2))
print(summary(lm_pos2)$coefficients)
print(summary(lm_pos2))
write(c("pos193",coefficients(lm_pos2)[2]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos193",summary(lm_pos2)$r.squared), sep=",", file="PTM_R2_Terms.csv", append=TRUE, ncol=2)

p_value_pos2 <- summary(lm_pos2)$coefficients[2,4]


lm_pos3 <- lm(DHC ~ pos258)

print("Position 258")
print(coefficients(lm_pos3))
print(summary(lm_pos3)$coefficients)
print(summary(lm_pos3))
write(c("pos258",coefficients(lm_pos3)[2]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos258",summary(lm_pos3)$r.squared), sep=",", file="PTM_R2_Terms.csv", append=TRUE, ncol=2)


p_value_pos3 <- summary(lm_pos3)$coefficients[2,4]


lm_pos4 <- lm(DHC ~ pos271)

print("Position 271")
print(coefficients(lm_pos4))
print(summary(lm_pos4)$coefficients)
print(summary(lm_pos4))
write(c("pos271",coefficients(lm_pos4)[2]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos271",summary(lm_pos4)$r.squared), sep=",", file="PTM_R2_Terms.csv", append=TRUE, ncol=2)

p_value_pos4 <- summary(lm_pos4)$coefficients[2,4]


lm_pos5 <- lm(DHC ~ pos273)

print("Position 273")
print(coefficients(lm_pos5))
print(summary(lm_pos5)$coefficients)
print(summary(lm_pos5))
write(c("pos273",coefficients(lm_pos5)[2]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos273",summary(lm_pos5)$r.squared), sep=",", file="PTM_R2_Terms.csv", append=TRUE, ncol=2)

p_value_pos5 <- summary(lm_pos5)$coefficients[2,4]


# Full First Order Model

lm_firstOrder <- lm(DHC ~ pos72 + pos193 + pos258 + pos271 + pos273)

print("Full First Order Model")
print(coefficients(lm_firstOrder))
print(summary(lm_firstOrder)$coefficients)
print(summary(lm_firstOrder))
write(c("First",summary(lm_firstOrder)$r.squared), sep=",", file="PTM_R2_Terms.csv", append=TRUE, ncol=2)
first_order_R2 <- summary(lm_firstOrder)$r.squared


p_value_firstOrder <- pf(summary(lm_firstOrder)$fstatistic[1], summary(lm_firstOrder)$fstatistic[2], summary(lm_firstOrder)$fstatistic[3], lower.tail = FALSE)


# Epistatic Interactions - Second Order

pos72pos193 <- pos72*pos193
pos72pos258 <- pos72*pos258
pos72pos271 <- pos72*pos271
pos72pos273 <- pos72*pos273
pos193pos258 <- pos193*pos258
pos193pos271 <- pos193*pos271
pos193pos273 <- pos193*pos273
pos258pos271 <- pos258*pos271
pos258pos273 <- pos258*pos273
pos271pos273 <- pos271*pos273

lm_secondOrder <- lm(DHC ~ pos72 + pos193 + pos258 + pos271 + pos273 + pos72pos193 + pos72pos258 + pos72pos271 + pos72pos273 + pos193pos258 + pos193pos271 + pos193pos273 + pos258pos271 + pos258pos273 + pos271pos273)

print("Full Second Order Model")
print(coefficients(lm_secondOrder))
print(summary(lm_secondOrder)$coefficients)
print(summary(lm_secondOrder))

write(c("pos72pos193",coefficients(lm_secondOrder)[7]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos72pos258",coefficients(lm_secondOrder)[8]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos72pos271",coefficients(lm_secondOrder)[9]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos72pos273",coefficients(lm_secondOrder)[10]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos193pos258",coefficients(lm_secondOrder)[11]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos193pos271",coefficients(lm_secondOrder)[12]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos193pos273",coefficients(lm_secondOrder)[13]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos258pos271",coefficients(lm_secondOrder)[14]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos258pos273",coefficients(lm_secondOrder)[15]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos271pos273",coefficients(lm_secondOrder)[16]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)


second_order_R2 <- summary(lm_secondOrder)$r.squared
Diff_second_order <- second_order_R2 - first_order_R2
write(c("Second",Diff_second_order), sep=",", file="PTM_R2_Terms.csv", append=TRUE, ncol=2)

lrtest_secondOrder <- lrtest(lm_secondOrder, lm_firstOrder)

print(lrtest_secondOrder)

SecondOrder_store <- summary(lrtest_secondOrder)

print(SecondOrder_store)

print(lrtest_secondOrder[2,5])

p_value_secondOrder <- lrtest_secondOrder[2,5]


# Epistatic Interactions - Third Order

pos72pos193pos258 <- pos72*pos193*pos258
pos72pos193pos271 <- pos72*pos193*pos271
pos72pos193pos273 <- pos72*pos193*pos273
pos72pos258pos271 <- pos72*pos258*pos271
pos72pos258pos273 <- pos72*pos258*pos273
pos72pos271pos273 <- pos72*pos271*pos273

pos193pos258pos271 <- pos193*pos258*pos271
pos193pos258pos273 <- pos193*pos258*pos273
pos193pos271pos273 <- pos193*pos271*pos273

pos258pos271pos273 <- pos258*pos271*pos273

lm_thirdOrder <- lm(DHC ~ pos72 + pos193 + pos258 + pos271 + pos273 + pos72pos193 + pos72pos258 + pos72pos271 + pos72pos273 + pos193pos258 + pos193pos271 + pos193pos273 + pos258pos271 + pos258pos273 + pos271pos273 + pos72pos193pos258 + pos72pos193pos271 + pos72pos193pos273 + pos72pos258pos271 + pos72pos258pos273 + pos72pos271pos273 + pos193pos258pos271 + pos193pos258pos273 + pos193pos271pos273 + pos258pos271pos273)

print("Full Third Order Model")
print(coefficients(lm_thirdOrder))
print(summary(lm_thirdOrder)$coefficients)
print(summary(lm_thirdOrder))

write(c("pos72pos193pos258",coefficients(lm_thirdOrder)[17]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos72pos193pos271",coefficients(lm_thirdOrder)[18]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos72pos193pos273",coefficients(lm_thirdOrder)[19]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos72pos258pos271",coefficients(lm_thirdOrder)[20]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos72pos258pos273",coefficients(lm_thirdOrder)[21]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos72pos271pos273",coefficients(lm_thirdOrder)[22]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos193pos258pos271",coefficients(lm_thirdOrder)[23]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos193pos258pos273",coefficients(lm_thirdOrder)[24]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos193pos271pos273",coefficients(lm_thirdOrder)[25]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos258pos271pos273",coefficients(lm_thirdOrder)[26]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)


third_order_R2 <- summary(lm_thirdOrder)$r.squared
Diff_third_order <- third_order_R2 - second_order_R2
write(c("Third",Diff_third_order), sep=",", file="PTM_R2_Terms.csv", append=TRUE, ncol=2)

lrtest_thirdOrder <- lrtest(lm_thirdOrder, lm_secondOrder)

print(lrtest_thirdOrder)

ThirdOrder_store <- summary(lrtest_thirdOrder)

p_value_thirdOrder <- lrtest_thirdOrder[2,5]


# Epistatic Interactions - Fourth Order

pos72pos193pos258pos271 <- pos72*pos193*pos258*pos271
pos72pos193pos258pos273 <- pos72*pos193*pos258*pos273
pos72pos193pos271pos273 <- pos72*pos193*pos271*pos273

pos72pos258pos271pos273 <- pos72*pos258*pos271*pos273

pos193pos258pos271pos273 <- pos193*pos258*pos271*pos273

lm_fourthOrder <- lm(DHC ~ pos72 + pos193 + pos258 + pos271 + pos273 + pos72pos193 + pos72pos258 + pos72pos271 + pos72pos273 + pos193pos258 + pos193pos271 + pos193pos273 + pos258pos271 + pos258pos273 + pos271pos273 + pos72pos193pos258 + pos72pos193pos271 + pos72pos193pos273 + pos72pos258pos271 + pos72pos258pos273 + pos72pos271pos273 + pos193pos258pos271 + pos193pos258pos273 + pos193pos271pos273 + pos258pos271pos273 + pos72pos193pos258pos271 + pos72pos193pos258pos273 + pos72pos193pos271pos273 + pos72pos258pos271pos273 + pos193pos258pos271pos273)

print("Full Fourth Order Model")
print(coefficients(lm_fourthOrder))
print(summary(lm_fourthOrder)$coefficients)
print(summary(lm_fourthOrder))

write(c("pos72pos193pos258pos271",coefficients(lm_fourthOrder)[27]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos72pos193pos258pos273",coefficients(lm_fourthOrder)[28]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos72pos193pos271pos273",coefficients(lm_fourthOrder)[29]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos72pos258pos271pos273",coefficients(lm_fourthOrder)[30]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)
write(c("pos193pos258pos271pos273",coefficients(lm_fourthOrder)[31]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)


fourth_order_R2 <- summary(lm_fourthOrder)$r.squared
Diff_fourth_order <- fourth_order_R2 - third_order_R2
write(c("Fourth",Diff_fourth_order), sep=",", file="PTM_R2_Terms.csv", append=TRUE, ncol=2)


lrtest_fourthOrder <- lrtest(lm_fourthOrder, lm_thirdOrder)

print(lrtest_fourthOrder)

FourthOrder_store <- summary(lrtest_fourthOrder)

p_value_fourthOrder <- lrtest_fourthOrder[2,5]


# Epistatic Interactions - Fifth Order

pos72pos193pos258pos271pos273 <- pos72*pos193*pos258*pos271*pos273

lm_fifthOrder <- lm(DHC ~ pos72 + pos193 + pos258 + pos271 + pos273 + pos72pos193 + pos72pos258 + pos72pos271 + pos72pos273 + pos193pos258 + pos193pos271 + pos193pos273 + pos258pos271 + pos258pos273 + pos271pos273 + pos72pos193pos258 + pos72pos193pos271 + pos72pos193pos273 + pos72pos258pos271 + pos72pos258pos273 + pos72pos271pos273 + pos193pos258pos271 + pos193pos258pos273 + pos193pos271pos273 + pos258pos271pos273 + pos72pos193pos258pos271 + pos72pos193pos258pos273 + pos72pos193pos271pos273 + pos72pos258pos271pos273 + pos193pos258pos271pos273 + pos72pos193pos258pos271pos273)

print("Full Fifth Order Model")
print(coefficients(lm_fifthOrder))
print(summary(lm_fifthOrder)$coefficients)
print(summary(lm_fifthOrder))

write(c("pos72pos193pos258pos271pos273",coefficients(lm_fifthOrder)[32]), sep=",", file="PTM_Terms.csv", append=TRUE, ncol=2)

fifth_order_R2 <- summary(lm_fifthOrder)$r.squared
Diff_fifth_order <- fifth_order_R2 - fourth_order_R2
write(c("Fifth",Diff_fifth_order), sep=",", file="PTM_R2_Terms.csv", append=TRUE, ncol=2)


lrtest_fifthOrder <- lrtest(lm_fifthOrder, lm_fourthOrder)

print(lrtest_fifthOrder)

FifthOrder_store <- summary(lrtest_fifthOrder)

p_value_fifthOrder <- lrtest_fifthOrder[2,5]

p_values_DHC <- c(p_value_pos1, p_value_pos2, p_value_pos3, p_value_pos4, p_value_pos5, p_value_firstOrder, p_value_secondOrder, p_value_thirdOrder, p_value_fourthOrder, p_value_fifthOrder)

print("P Values for DHC")
print(cbind(c("Pos1", "Pos2", "Pos3", "Pos4", "Pos5", "First Order", "Second Order", "Third Order", "Fourth Order", "Fifth Order"), p_values_DHC))