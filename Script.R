
# Script to do all analyses in Simpson et al.
# Written by AM Senior at the University of Sydney, 2024
# Chunks of Code are broadly split up by Figure and associated analyses, each chunk is independent of the others and can be run on its own

# Files required in the working directory are
# 1) Hu_et_al_Intakes.csv - the intake in grams and bwt of animals in their study
# 2) Hu_et_al_Diets.csv - the diet composoitions for that study - intakes of specific macros etc. are derived from these values and intake in grams
# 3) Solon_Biet_et_al_Intakes.csv - the intake in grams and bwt of animals in their study
# 4) Solon_Biet_et_al_Diets.csv - the diet composoitions for that study - intakes of specific macros etc. are derived from these values and intake in grams
# 5) Solon_Biet_et_al_Intakes6mo.csv - the intake in grams and bwt of animals in their study for the first 6 months
# 6 - 9) BALB-c.csv, C3H.csv, DBA2.csv, FVB.csv - intakes for additional mouse strains tested by Hu et al.
# 10) Hu_et_al_Diets_in_S1.csv - the diet composoitions for that study as stated in thier supplementary materials

############################################################
######################## Figure 1 ##########################
############################################################

# Clear any old objects
rm(list=ls())

# Load the libraries
library(compositions)

####

# size of points
psize<-1.65

pdf("Figure1.pdf", height=10, width=10)

par(mfrow=c(2,2), mar=c(5,5,5,5), cex.lab=1.6, cex.axis=1.4)

# Create EMT for Solon=Biet et al.
data<-read.csv("Solon_Biet_et_al_Diets.csv")
data$P<-data$diet.P / (data$diet.P + data$diet.C + data$diet.F) * 100
data$C<-data$diet.C / (data$diet.P + data$diet.C + data$diet.F) * 100
data$F<-data$diet.F / (data$diet.P + data$diet.C + data$diet.F) * 100

data.aitch<-acomp(data, parts=c("F", "P", "C"))
plot.acomp(data.aitch, col=1, cex=psize, labels=c("Fat", "Protein", "Carbohydrate"), pch=16)
isoPortionLines()

text(53, 53, "%Fat", srt=-45, cex=1.6)
mtext("A", at=-0.1, line=1, cex=2, font=2)
mtext("Solon-Biet et al. 2014", line=0, cex=1.25)


# Create EMT for Hu et al.
data<-read.csv("Hu_et_al_Diets.csv")
data$P<-data$percent.P / (data$percent.P + data$percent.C + data$percent.F) * 100
data$C<-data$percent.C / (data$percent.P + data$percent.C + data$percent.F) * 100
data$F<-data$percent.F / (data$percent.P + data$percent.C + data$percent.F) * 100

data.aitch<-acomp(data, parts=c("F", "P", "C"))
plot.acomp(data.aitch, col=1, cex=psize, labels=c("Fat", "Protein", "Carbohydrate"), pch=16)
isoPortionLines()

text(53, 53, "%Fat", srt=-45, cex=1.6)
mtext("B", at=-0.1, line=1, cex=2, font=2)
mtext("Hu et al. 2018", line=0, cex=1.25)


# Fat against energy density in Hu et al.
plot(data$percent.F, data$energy.density * 4.184, xlab="% Fat", ylab="Energy Density (kJ/g)", pch=16, cex=psize)
mtext("C", at=1, line=2, cex=2, font=2)

# Fat and carbohydrate against energy density in Hu et al.
plot(data$percent.F + data$percent.C, data$energy.density * 4.184, xlab="% Non-Protein Energy", ylab="Energy Density (kJ/g)", pch=16, cex=psize)
mtext("D", at=67, line=2, cex=2, font=2)


dev.off()

############################################################
######################## Figure 2 ##########################
############################################################

# Clear any old objects
rm(list=ls())

# Load the libraries
library(mgcv)
library(sp)
library(fields)

# A function created to find the outer perimeter over which the surface should be fitted
findConvex<-function(x,y,rgnames,res=101){
	hull<-cbind(x,y)[chull(cbind(x,y)),]
	px<-pretty(x)
	py<-pretty(y)
	x.new<-seq(min(px),max(px),len=res)
	y.new<-seq(min(py),max(py),len=res)
	ingrid<-as.data.frame(expand.grid(x.new,y.new))                                                              
	Fgrid<-ingrid
	Fgrid[(point.in.polygon(ingrid[,1], ingrid[,2], hull[,1],hull[,2])==0),]<-NA
	names(Fgrid)<-rgnames
	return(Fgrid)
}

# For this analysis we need to know the kJ of P and kJ of non P per gram in each animals diet, as well as energy intake

####

# Data from Hu et al. paper
# Read the animal data
data<-read.csv("Hu_et_al_Intakes.csv")
head(data)

# Read in the dietary data
diets<-read.csv("Hu_et_al_Diets.csv")

# Convert the calories per gram to kJs
diets$kJ.g<-diets$energy.density * 4.184

# match up the diets to the intakes and find the compositon and energy density of each animals intake
data$kJ.g<-diets$kJ.g[match(data$diet, diets$diet)]
data$per.P<-diets$percent.P[match(data$diet, diets$diet)]
data$per.C<-diets$percent.C[match(data$diet, diets$diet)]
data$per.F<-diets$percent.F[match(data$diet, diets$diet)]

# Calculate the amount of P, C and F in the diets
data$diet.P<-data$kJ.g * (data$per.P / 100)
data$diet.C<-data$kJ.g * (data$per.C / 100)
data$diet.F<-data$kJ.g * (data$per.F / 100)
data$diet.nonP<-data$diet.C + data$diet.F

# Calculate energy intakes
data$intake.E<-data$intake.gram * data$kJ.g

# Save the dataset in to a list and remove it
datasets<-list()
datasets[[1]]<-data
rm(data)
rm(diets)

# Load the data from the GF study
data<-read.csv("Solon_Biet_et_al_Intakes.csv")
head(data)

# Read in the dietary data
diets<-read.csv("Solon_Biet_et_al_Diets.csv")

# match up the diets to the intakes and find the Kjs of each macro per gram of each animals diet
data$diet.P<-diets$diet.P[match(data$diet, diets$diet)]
data$diet.C<-diets$diet.C[match(data$diet, diets$diet)]
data$diet.F<-diets$diet.F[match(data$diet, diets$diet)]
data$diet.nonP<-data$diet.C + data$diet.F

# Energy content per gram
data$kJ.g<-data$diet.P + data$diet.C + data$diet.F

# Calculate energy intakes
data$intake.E<-data$intake.gram * data$kJ.g

# Save the dataset and remove old object
datasets[[2]]<-data
rm(data)

# Plotting

# Create plot for intakes
pdf("Figure2.pdf", height=11, width=10)

par(mar=c(5,5,6,1), cex.axis=1.3)
layout(rbind(c(1,3), c(2,4)))

# Plot surfaces from Hu et al data alone

# Get the right dataset
data<-datasets[[1]]

# Set the resolution of the surface
surface.resolution<-501

# How many values to round surface
round.surf<-3

# This specifies the color scheme for surface - it is actually a function that returns a function
rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")

# How many different colours should we use on the plot
no.cols<-256

# Get the colors to use from the pallette specified above
map<-rgb.palette(no.cols)

# How many levels should there be on the surface
nlev<-4

# surfaces in order of presentation for x, y, and median slices
surfaces.list<-c("diet.P", "diet.nonP")

# Labels for each panel
labels<-c("Protein (kJ/g)", "Non Protein (kJ/g)")

# Specify the model formula to use
model.form<-as.formula("outcome.j ~ s(diet.P, diet.nonP, k=k)")

# Restrict k
k<-12

# List the outcomes to model
outcomes<-c("intake.gram", "intake.E")
names.outcomes<-c("Food Intake (g)", "Energy Intake (kJ)")

# Titles for the different outcomes
titles<-c("Hu et al. 2018", "Solon-Biet et al. 2014")

# Letters for the panels
letters<-c("A", "C")

# list to hold the models
models<-list()

# Open the loop for the outcomes
for(j in 1:length(outcomes)){
	
	# find thr jth outcome
	data$outcome.j<-data[,outcomes[j]]
	
	# Fit the GAM
	model<-gam(model.form, data=data, family=scat())
	
	# Predicted values for surface
	Pred.Values<-findConvex(data[,surfaces.list[1]], data[,surfaces.list[2]], c(surfaces.list[1], surfaces.list[2]), res=surface.resolution)
		
	out<-predict(model, newdata = Pred.Values, type = "response")
	surf<-matrix(out, nrow=sqrt(dim(Pred.Values)[1]))
	surf<-round(surf, round.surf)
		
	# Pretty comes up with nice values of x and y over which to fit
	px<-pretty(data[,surfaces.list[1]])
	py<-pretty(data[,surfaces.list[2]])
		
	# Uses px and py to generate the x and y axes
	x.new<-seq(min(px), max(px), len=surface.resolution)
	y.new<-seq(min(py), max(py), len=surface.resolution)
		
	# We need to know the minimum and maximum predicted values to make the scale of the plots sensible
	mn<-min(out, na.rm=T)
	mx<-max(out, na.rm=T)
	locs<-(range(out, na.rm=TRUE) - mn) / (mx-mn) * no.cols
				
	# Actually plots the surface using all of the above info above
	image(x.new, y.new, surf, col=map[locs[1]:locs[2]], xlab="", ylab="", axes=FALSE, main="", xlim=c(0, 12), ylim=c(0, 25))
	# Adds some axes
	axis(1)
	axis(2)
	# Adds a contour over the top (can add a title using main)
	contour(x.new, y.new, surf, add=TRUE, levels=pretty(range(mn,mx), nlev), labcex=0.8)
		
	# Add the labels	
	mtext(labels[1], side=1, cex=1.25, line=3)
	mtext(labels[2], side=2, cex=1.25, line=3)		
	mtext(titles[1], line=-1, cex=1.5)
	mtext(letters[j], font=2, at=-1, line=1.5, cex=2.5)
	mtext(names.outcomes[j], font=2, at=13, line=4, cex=1.75)
	
	# Save the model and remove it
	models[[j]]<-model
	rm(model)

}

# Now do exactly the same as above but with the GF data
# Get the right dataset
data<-datasets[[2]]

# Letters for the panels
letters<-c("B", "D")

# Open the loop for the outcomes
for(j in 1:length(outcomes)){
	
	# find thr jth outcome
	data$outcome.j<-data[,outcomes[j]]
	
	# Fit the GAM
	model<-gam(model.form, data=data, family=scat())
	
	# Predicted values for surface
	Pred.Values<-findConvex(data[,surfaces.list[1]], data[,surfaces.list[2]], c(surfaces.list[1], surfaces.list[2]), res=surface.resolution)
		
	out<-predict(model, newdata = Pred.Values, type = "response")
	surf<-matrix(out, nrow=sqrt(dim(Pred.Values)[1]))
	surf<-round(surf, round.surf)
		
	# Pretty comes up with nice values of x and y over which to fit
	px<-pretty(data[,surfaces.list[1]])
	py<-pretty(data[,surfaces.list[2]])
		
	# Uses px and py to generate the x and y axes
	x.new<-seq(min(px), max(px), len=surface.resolution)
	y.new<-seq(min(py), max(py), len=surface.resolution)
		
	# We need to know the minimum and maximum predicted values to make the scale of the plots sensible
	mn<-min(out, na.rm=T)
	mx<-max(out, na.rm=T)
	locs<-(range(out, na.rm=TRUE) - mn) / (mx-mn) * no.cols
				
	# Actually plots the surface using all of the above info above
	image(x.new, y.new, surf, col=map[locs[1]:locs[2]], xlab="", ylab="", axes=FALSE, main="", xlim=c(0, 12), ylim=c(0, 25))
	# Adds some axes
	axis(1)
	axis(2)
	# Adds a contour over the top (can add a title using main)
	contour(x.new, y.new, surf, add=TRUE, levels=pretty(range(mn,mx), nlev), labcex=0.8)
		
	# Add the labels	
	mtext(labels[1], side=1, cex=1.25, line=3)
	mtext(labels[2], side=2, cex=1.25, line=3)		
	mtext(titles[2], line=-1, cex=1.5)
	mtext(letters[j], font=2, at=-1, line=1.5, cex=2.5)
	
	# Save the model and remove it
	models[[2+j]]<-model
	rm(model)
}

dev.off()

# Table to hold results
results<-as.data.frame(array(NA, c(4,7)))
names(results)<-c("Outcome", "Study", "Coef.", "edf", "rdf", "F", "p")
results$Outcome<-c(rep("Mass Intake", 2), rep("Energy Intake", 2))
results$Study<-rep(c(rep("Hu et al.", 1), rep("Solon-Biet et al.", 1)), 2)

# Order of models in table
table.order<-c(1,3,2,4)
counter<-1

# Add in the model results
for(i in 1:4){
	
	# Pull out the model and save in the table then remove it
	model<-models[[table.order[i]]]
	results[i,c(4:7)]<-summary(model)$s.table
	results[i,3]<-rownames(summary(model)$s.table)
	rm(model)
	
}

# Round the figures and write a results table 
results[,c(4:6)]<-round(results[,c(4:6)], 3)
write.table(results, file="TableS3.pt1.csv", sep=",", row.names=F, col.names=names(results))

############################################################
######################## Figure S2 #########################
############################################################

# Clear any old objects
rm(list=ls())

# Load the libraries
library(mgcv)
library(sp)
library(fields)

# A function created to find the outer perimeter over which the surface should be fitted
findConvex<-function(x,y,rgnames,res=101){
	hull<-cbind(x,y)[chull(cbind(x,y)),]
	px<-pretty(x)
	py<-pretty(y)
	x.new<-seq(min(px),max(px),len=res)
	y.new<-seq(min(py),max(py),len=res)
	ingrid<-as.data.frame(expand.grid(x.new,y.new))                                                              
	Fgrid<-ingrid
	Fgrid[(point.in.polygon(ingrid[,1], ingrid[,2], hull[,1],hull[,2])==0),]<-NA
	names(Fgrid)<-rgnames
	return(Fgrid)
}

# For this analysis we need to calculate the intakes of energy for each animal in each strain, as well as the P nonPkJ/gram of diet for that aniaml

# Data from Hu et al.

strains<-c("BALB-c", "C3H", "DBA2", "FVB")

# Read the animal data
data<-read.csv(paste(strains[1], ".csv", sep=""))
data$strain<-strains[1]
for(i in 2:4){
	data.i<-read.csv(paste(strains[i], ".csv", sep=""))
	data.i$strain<-strains[i]
	data<-rbind(data, data.i)	
}

# Read in C57bl6 data and trim down to comparable diets only
c57<-read.csv("Hu_et_al_Intakes.csv")
c57<-c57[-which(is.na(match(c57$diet, data$Diet.Group)) == T),]
c57$strain<-"C57Bl6"

# Chuck out a few columns we do not need
c57<-c57[,c(2,4,5)]

# Add to the others
names(data)<-names(c57)
data<-rbind(data, c57)

# Drop a few rows with unkown diets - weird
data<-data[-which(is.na(data$diet) == T),]

# Read in the dietary data
diets<-read.csv("Hu_et_al_Diets.csv")

# Convert the calories per gram to kJs
diets$kJ.g<-diets$energy.density * 4.184

# match up the diets to the intakes and find the compositon and energy density of each animals intake
data$kJ.g<-diets$kJ.g[match(data$diet, diets$diet)]
data$per.P<-diets$percent.P[match(data$diet, diets$diet)]
data$per.C<-diets$percent.C[match(data$diet, diets$diet)]
data$per.F<-diets$percent.F[match(data$diet, diets$diet)]

# Calculate energy intakes
data$intake.E<-data$intake.g * data$kJ.g

# Calculate the amount of P, C and F in the diets
data$diet.P<-data$kJ.g * (data$per.P / 100)
data$diet.C<-data$kJ.g * (data$per.C / 100)
data$diet.F<-data$kJ.g * (data$per.F / 100)
data$diet.nonP<-data$diet.C + data$diet.F

# Create plot for intakes
pdf("FigureS2.pdf", height=15, width=18)
layout(cbind(seq(1, 6, 2), seq(7, 12, 2), seq(2, 6, 2), seq(8, 12, 2)))
par(mar=c(5,5,5,5), cex.axis=1.5)

# Set the resolution of the surface
surface.resolution<-501

# How many values to round surface
round.surf<-3

# This specifies the color scheme for surface - it is actually a function that returns a function
rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")

# How many different colours should we use on the plot
no.cols<-256

# Get the colors to use from the pallette specified above
map<-rgb.palette(no.cols)

# How many levels should there be on the surface
nlev<-4

# surfaces in order of presentation for x, y, and median slices
surfaces.list<-c("diet.P", "diet.nonP")

# Labels for each panel
labels<-c("Protein (kJ/g)", "Non Protein (kJ/g)")

# Specify the model formula to use
model.form<-as.formula("outcome.j ~ s(diet.P, diet.nonP, k=k)")

# Restrict k
k<-12

# List the outcomes to model
outcomes<-c("intake.gram", "intake.E")

# Titles for the different outcomes
titles<-c("Food Intake (g)", "Energy Intake (kJ)")

# List for models
models<-list()

# Start by doing the pooled data

# Open the loop for the outcomes
for(j in 1:length(outcomes)){
	
	# find thr jth outcome
	data$outcome.j<-data[,outcomes[j]]
	
	# Fit the GAM
	model<-gam(model.form, data=data, family=scat())
	
	# Predicted values for surface
	Pred.Values<-findConvex(data[,surfaces.list[1]], data[,surfaces.list[2]], c(surfaces.list[1], surfaces.list[2]), res=surface.resolution)
		
	out<-predict(model, newdata = Pred.Values, type = "response")
	surf<-matrix(out, nrow=sqrt(dim(Pred.Values)[1]))
	surf<-round(surf, round.surf)
		
	# Pretty comes up with nice values of x and y over which to fit
	px<-pretty(data[,surfaces.list[1]])
	py<-pretty(data[,surfaces.list[2]])
		
	# Uses px and py to generate the x and y axes
	x.new<-seq(min(px), max(px), len=surface.resolution)
	y.new<-seq(min(py), max(py), len=surface.resolution)
		
	# We need to know the minimum and maximum predicted values to make the scale of the plots sensible
	mn<-min(out, na.rm=T)
	mx<-max(out, na.rm=T)
	locs<-(range(out, na.rm=TRUE) - mn) / (mx-mn) * no.cols
				
	# Actually plots the surface using all of the above info above
	image(x.new, y.new, surf, col=map[locs[1]:locs[2]], xlab="", ylab="", axes=FALSE, xlim=c(0, 8), ylim=c(12, 24))
	# Adds some axes
	axis(1)
	axis(2)
	# Adds a contour over the top (can add a title using main)
	contour(x.new, y.new, surf, add=TRUE, levels=pretty(range(mn,mx), nlev), labcex=0.8)
		
	# Add the labels	
	mtext(labels[1], side=1, cex=1.4, line=3)
	mtext(labels[2], side=2, cex=1.4, line=3)		
	mtext("All Strains", line=-1, cex=1.7)
	mtext(titles[j], at=10, line=2, cex=1.8, font=2)
	
	# Save the model
	models[[j]]<-model
	names(models)[[j]]<-paste(outcomes[j], "All Strains", sep=" - ")
}


# Repeat above for individual strains
strains<-c("C57Bl6", "BALB-c", "C3H", "DBA2", "FVB")
counter<-3

for(strain in 1:5){
	
	data.strain<-data[which(data$strain == strains[strain]),]
	
	## Plot surfaces from speakman data alone
	
	# Set the resolution of the surface
	surface.resolution<-501
	
	# How many values to round surface
	round.surf<-3
	
	# This specifies the color scheme for surface - it is actually a function that returns a function
	rgb.palette<-colorRampPalette(c("blue","cyan","yellow","red"), space="Lab", interpolate="linear")
	
	# How many different colours should we use on the plot
	no.cols<-256
	
	# Get the colors to use from the pallette specified above
	map<-rgb.palette(no.cols)
	
	# How many levels should there be on the surface
	nlev<-4
	
	# surfaces in order of presentation for x, y, and median slices
	surfaces.list<-c("diet.P", "diet.nonP")
	
	# Labels for each panel
	labels<-c("Protein (kJ/g)", "Non Protein (kJ/g)")
	
	# Specify the model formula to use
	model.form<-as.formula("outcome.j ~ s(diet.P, diet.nonP, k=k)")
	
	# Restrict k
	k<-12
	
	# List the outcomes to model
	outcomes<-c("intake.gram", "intake.E")
	
	# Open the loop for the outcomes
	for(j in 1:length(outcomes)){
		
		# find thr jth outcome
		data.strain$outcome.j<-data.strain[,outcomes[j]]
		
		# Fit the GAM
		model<-gam(model.form, data=data.strain, family=scat())
		
		# Predicted values for surface
		Pred.Values<-findConvex(data.strain[,surfaces.list[1]], data.strain[,surfaces.list[2]], c(surfaces.list[1], surfaces.list[2]), res=surface.resolution)
			
		out<-predict(model, newdata = Pred.Values, type = "response")
		surf<-matrix(out, nrow=sqrt(dim(Pred.Values)[1]))
		surf<-round(surf, round.surf)
			
		# Pretty comes up with nice values of x and y over which to fit
		px<-pretty(data[,surfaces.list[1]])
		py<-pretty(data[,surfaces.list[2]])
			
		# Uses px and py to generate the x and y axes
		x.new<-seq(min(px), max(px), len=surface.resolution)
		y.new<-seq(min(py), max(py), len=surface.resolution)
			
		# We need to know the minimum and maximum predicted values to make the scale of the plots sensible
		mn<-min(out, na.rm=T)
		mx<-max(out, na.rm=T)
		locs<-(range(out, na.rm=TRUE) - mn) / (mx-mn) * no.cols
					
		# Actually plots the surface using all of the above info above
		image(x.new, y.new, surf, col=map[locs[1]:locs[2]], xlab="", ylab="", axes=FALSE, xlim=c(0, 8), ylim=c(12, 24))
		# Adds some axes
		axis(1)
		axis(2)
		# Adds a contour over the top (can add a title using main)
		contour(x.new, y.new, surf, add=TRUE, levels=pretty(range(mn,mx), nlev), labcex=0.8)
			
		# Add the labels	
		mtext(labels[1], side=1, cex=1.4, line=3)
		mtext(labels[2], side=2, cex=1.4, line=3)		
		mtext(strains[strain], line=-1, cex=1.7)
		
		# Save the model
		models[[counter]]<-model
		names(models)[[counter]]<-paste(outcomes[j], strains[strain], sep=" - ")
		counter<-counter+1
	}

}

dev.off()


# Table to hold results
results<-as.data.frame(array(NA, c(2*6, 7)))
names(results)<-c("Outcome", "Strain", "Coef.", "edf", "rdf", "F", "p")
results$Outcome<-c(rep("Mass Intake", 6), rep("Energy Intake", 6))
results$Strain<-rep(c("All Strains", "C57Bl6", "BALB-c", "C3H", "DBA2", "FVB"), 2)

# Order of models in table
table.order<-c(seq(1, length(models), 2), seq(2, length(models), 2))
counter<-1

# Add in the model results
for(i in 1:12){
	
	# Pull out the model and save in the table then remove it
	model<-models[[table.order[i]]]
	results[i ,c(4:7)]<-summary(model)$s.table
	results[i ,3]<-rownames(summary(model)$s.table)
	rm(model)
	
}

# Round the figures and write the table
results[,c(4:6)]<-round(results[,c(4:6)], 3)
write.table(results, file="TableS3.pt2.csv", sep=",", row.names=F, col.names=names(results))

############################################################
######################## Figure S3 #########################
############################################################

# Clear any old objects
rm(list=ls())

# For this analysis we need to calculate energy intake on of each animal on each diet

####

# Data from Hu et al paper
# Read the animal data
data<-read.csv("Hu_et_al_Intakes.csv")
head(data)

# Read in the dietary data
diets<-read.csv("Hu_et_al_Diets.csv")

# Convert the calories per gram to kJs
diets$kJ.g<-diets$energy.density * 4.184

# match up the diets to the intakes and find the compositon and energy density of each animals intake
data$kJ.g<-diets$kJ.g[match(data$diet, diets$diet)]

# Calculate energy intakes
data$intake.E<-data$intake.gram * data$kJ.g

# Add to the list
datasets<-list()
datasets[[1]]<-data

# Do the same for Solon-Biet
# Read the animal data
data<-read.csv("Solon_Biet_et_al_Intakes.csv")
head(data)

# Read in the dietary data
diets<-read.csv("Solon_Biet_et_al_Diets.csv")

# match up the diets to the intakes and find the Kjs of each macro per gram of each animals diet
data$diet.P<-diets$diet.P[match(data$diet, diets$diet)]
data$diet.C<-diets$diet.C[match(data$diet, diets$diet)]
data$diet.F<-diets$diet.F[match(data$diet, diets$diet)]

# Calculate the kJs per gram
data$kJ.g<-data$diet.P + data$diet.C + data$diet.F

# Calculate energy intakes
data$intake.E<-data$intake.gram * data$kJ.g

# Add to the data list and remove
datasets[[2]]<-data
rm(data)

# Repeat also for 6Mo data in Solon-Biet
data<-read.csv("Solon_Biet_et_al_Intakes_6mo.csv")
head(data)

# match up the diets to the intakes and find the Kjs of each macro per gram of each animals diet
data$diet.P<-diets$diet.P[match(data$diet, diets$diet)]
data$diet.C<-diets$diet.C[match(data$diet, diets$diet)]
data$diet.F<-diets$diet.F[match(data$diet, diets$diet)]

# Calculate the kJs per gram
data$kJ.g<-data$diet.P + data$diet.C + data$diet.F

# Calculate energy intakes
data$intake.E<-data$intake.gram * data$kJ.g

# Add to the data list and remove
datasets[[3]]<-data
rm(data)

# Plotting 

# Always use the Hu et al. data to plot first
data<-datasets[[1]]

# Create plot for intakes
pdf("FigureS3.pdf", height=5, width=15)

pcar<-16
psize<-0.75

par(mfrow=c(1,3), mar=c(5,5,5,5), cex.lab=1.6, cex.axis=1.4)

# Use the full gf dataset
gf_data<-datasets[[2]]
plot(data$intake.E, as.numeric(data$b.wt), ylim=c(0, 60), xlim=c(0, 110), ylab="Body Mass (g)", xlab="Energy Intake (kJ)", pch=pcar, cex=psize)
points(gf_data$intake.E, gf_data$bwt, col=2, pch=pcar, cex=psize)
mtext("A", at=-5, font=2, cex=2, line=2)
legend(0, 60, c("Hu et al. 2018", "Solon-Biet et al. 2014"), pch=16, col=c(1,2), cex=1.4)

# Trim the gf data to males
gf_data<-gf_data[which(gf_data$sex == "M"),]
plot(data$intake.E, as.numeric(data$b.wt), ylim=c(0, 60), xlim=c(0, 110), ylab="Body Mass (g)", xlab="Energy Intake (kJ)", main="", pch=pcar, cex=psize)
points(gf_data$intake.E, gf_data$bwt, col=2, pch=pcar, cex=psize)
mtext("B", at=-5, font=2, cex=2, line=2)

# Use the data from the GF study at 6 months
gf_data<-datasets[[3]]
plot(data$intake.E, as.numeric(data$b.wt), ylim=c(0, 60), xlim=c(0, 110), ylab="Body Mass (g)", xlab="Energy Intake (kJ)", main="", pch=pcar, cex=psize)
points(gf_data$intake.E, gf_data$bwt, col=2, pch=pcar, cex=psize)
mtext("C", at=-5, font=2, cex=2, line=2)

dev.off()

############################################################
######################## Figure 3 ##########################
############################################################

# Clear any old objects
rm(list=ls())

# For this analysis we want to compare % energy from protein in the diet and intake in grams and energy
# And in only those diets with approximately equivelent energetic value

####

# Data from Hu et al. paper
# Read the animal data
data<-read.csv("Hu_et_al_Intakes.csv")
head(data)

# Read in the dietary data
diets<-read.csv("Hu_et_al_Diets.csv")

# Convert the calories per gram to kJs
diets$kJ.g<-diets$energy.density * 4.184

# match up the diets to the intakes and find the compositon and energy density of each animals intake
data$kJ.g<-diets$kJ.g[match(data$diet, diets$diet)]
data$per.P<-diets$percent.P[match(data$diet, diets$diet)]

# Constrain dataset to those less than 18kJ/g
data<-data[which(data$kJ.g < 18),]

# Calculate energy intakes
data$intake.E<-data$intake.gram * data$kJ.g

# Save the data to a list and remove it
datasets<-list()
datasets[[1]]<-data
rm(data)
rm(diets)

# Repeat for gf stusy
data<-read.csv("Solon_Biet_et_al_Intakes.csv")
head(data)

# Read in the dietary data
diets<-read.csv("Solon_Biet_et_al_Diets.csv")

# match up the diets to the intakes and find the Kjs of each macro per gram of each animals diet
data$diet.P<-diets$diet.P[match(data$diet, diets$diet)]
data$diet.C<-diets$diet.C[match(data$diet, diets$diet)]
data$diet.F<-diets$diet.F[match(data$diet, diets$diet)]
data$kJ.g<-data$diet.P + data$diet.C + data$diet.F

# Calculate percentage dietary protein
data$per.P<-data$diet.P / (data$diet.P + data$diet.F + data$diet.C) * 100

# Constrain dataset to those with high energy
data<-data[which(data$energy == "HIGH"),]

# Calculate energy intakes
data$intake.E<-data$intake.gram * data$kJ.g

# Save the dataset and remove
datasets[[2]]<-data
rm(data)
rm(diets)

# Plotting

# Create plot for intakes
pdf("Figure3.pdf", height=10, width=10)
par(mar=c(5,5,5,4), cex.lab=1.5, cex.axis=1.3)
layout(rbind(c(1,3), c(2,4)))

# Get the right dataset
data<-datasets[[1]]

# Labels for each panel
labels<-c("Protein (kJ/g)", "Non Protein (kJ/g)")

# List the outcomes to model
outcomes<-c("intake.gram", "intake.E")

# Titles for the different outcomes
titles<-c("Hu et al. 2018", "Solon-Biet et al. 2014")

# Y labels and limits
ylabs<-c("Food Intake (g)", "Energy Intake (kJ)")
ylims<-list()
ylims[[1]]<-c(0, 6)
ylims[[2]]<-c(0, 100)

# Letters for the panels
letters<-c("A", "C")

# list to hold the models
models<-list()

# Open the loop for the outcomes
for(j in 1:length(outcomes)){
	
	# find thr jth outcome
	data$outcome.j<-data[,outcomes[j]]
	
	# Fit the LM
	model<-lm(outcome.j ~ per.P, data=data)
	
	# Plot the data
	plot(data$per.P, data$outcome.j, xlab="% Protein", ylab=ylabs[j], xlim=c(0, 100), pch=16, cex=0.75, ylim=ylims[[j]])
	
	# Add the predictions - use a dashed line where there is no data
	pred.x<-seq(min(data$per.P), max(data$per.P), 0.1)
	y<-model$coef[1] + model$coef[2] * pred.x
	lines(pred.x, y)
	pred.x<-seq(0, min(data$per.P), 0.1)
	y<-model$coef[1] + model$coef[2] * pred.x
	lines(pred.x, y, lty=2)	
	pred.x<-seq(max(data$per.P), 100, 0.1)
	y<-model$coef[1] + model$coef[2] * pred.x
	lines(pred.x, y, lty=2)	
		
	# Add the labels	
	mtext(titles[1], line=-2, cex=1.5)
	mtext(letters[j], font=2, at=-1, line=1.5, cex=2.5)
	
	# Save the model and remove it
	models[[j]]<-model
	rm(model)

}

# remove the old data
rm(data)

# Now do exactly the same as above but with the GF data
# Get the right dataset
data<-datasets[[2]]

letters<-c("B", "D")

# Open the loop for the outcomes
for(j in 1:length(outcomes)){
	
	# find thr jth outcome
	data$outcome.j<-data[,outcomes[j]]
	
	# Fit the LM
	model<-lm(outcome.j ~ per.P, data=data)
	
	# Plot the data
	plot(data$per.P, data$outcome.j, xlab="% Protein", ylab=ylabs[j], xlim=c(0, 100), pch=16, cex=0.75, ylim=ylims[[j]])
	
	# Add the predictions - use a dashed line where there is no data
	pred.x<-seq(min(data$per.P), max(data$per.P), 0.1)
	y<-model$coef[1] + model$coef[2] * pred.x
	lines(pred.x, y)
	pred.x<-seq(0, min(data$per.P), 0.1)
	y<-model$coef[1] + model$coef[2] * pred.x
	lines(pred.x, y, lty=2)	
	pred.x<-seq(max(data$per.P), 100, 0.1)
	y<-model$coef[1] + model$coef[2] * pred.x
	lines(pred.x, y, lty=2)	
			
	# Add the labels	
	mtext(titles[2], line=-2, cex=1.5)
	mtext(letters[j], font=2, at=-1, line=1.5, cex=2.5)
	
	# Save the model and remove it
	models[[j+2]]<-model
	rm(model)

}

dev.off()

# Table to hold results
results<-as.data.frame(array(NA, c(8,7)))
names(results)<-c("Outcome", "Study", "Coef.", "Est.", "SE", "t", "p")
results$Outcome<-c(rep("Mass Intake", 4), rep("Energy Intake", 4))
results$Study<-rep(c(rep("Hu et al.", 2), rep("Solon-Biet et al.", 2)), 2)
results$Coef.<-rep(c("Intercept", "% Protein"), 4)

# Order of models in table
table.order<-c(1,3,2,4)
counter<-1

# Add in the model results
for(i in 1:4){
	
	# Pull out the model and save in the table then remove it
	model<-models[[table.order[i]]]
	results[c(counter:(counter+1)),c(4:7)]<-summary(model)$coefficients
	counter<-counter+2	
	rm(model)
	
}

# Round the figures and write the results
results[,c(4:6)]<-round(results[,c(4:6)], 3)
write.table(results, file="TableS4.csv", sep=",", row.names=F, col.names=names(results))


# For the model to contrast the slope of the two studies
pooled.data<-datasets[[1]][,c("intake.gram", "intake.E", "per.P")]
pooled.data$study<-as.factor("Hu_et_al")
data<-datasets[[2]][,c("intake.gram", "intake.E", "per.P")]
data$study<-as.factor("Solon_Biet_et_al")
pooled.data<-rbind(pooled.data, data)

# LM food intake in grams
model<-lm(intake.gram ~ per.P * study, data=pooled.data)
summary(model)

# LM food intake in energy
model<-lm(intake.E ~ per.P * study, data=pooled.data)
summary(model)

############################################################
######################## Figure 4 ##########################
############################################################

rm(list=ls())

# For this analysis we need to the proportion of protein per gram by dry weight in each diet

####

# Data from Hu et al. paper
# Read the animal data
data<-read.csv("Hu_et_al_Intakes.csv")
head(data)

# Read in the dietary data
diets<-read.csv("Hu_et_al_Diets.csv")

# Convert the calories per gram to kJs
diets$kJ.g<-diets$energy.density * 4.184

# match up the diets to the intakes and find the compositon and energy density of each animals intake
data$kJ.g<-diets$kJ.g[match(data$diet, diets$diet)]

# Work out how much protein is in the diet by grams
data$diet.g.P<-diets$Dry.P[match(data$diet, diets$diet)]

### Save the dataset in a list
datasets<-list()
datasets[[1]]<-data
rm(data)
rm(diets)

# Load the data from the GF study
# Read the animal data
data<-read.csv("Solon_Biet_et_al_Intakes.csv")
head(data)

# Read in the dietary data
diets<-read.csv("Solon_Biet_et_al_Diets.csv")

# Match up the dry weight by gram from the diet data
data$diet.g.P<-diets$Dry.P[match(data$diet, diets$diet)]

# match up the diets to the intakes and find the Kjs of each macro per gram of each animals diet
data$diet.P<-diets$diet.P[match(data$diet, diets$diet)]
data$diet.C<-diets$diet.C[match(data$diet, diets$diet)]
data$diet.F<-diets$diet.F[match(data$diet, diets$diet)]
data$kJ.g<-data$diet.P + data$diet.C + data$diet.F

# Save the data and remove objects
datasets[[2]]<-data
rm(data)
rm(diets)

# Plotting

pdf("Figure4.pdf", width=10, height=10)

# Aesthetics
pcars<-16
psize<-0.75
par(mfrow=c(2,2), mar=c(5,5,5,4), cex.lab=1.3, cex.axis=1.2)

# Titles for the figures
titles<-c("Hu et al. 2018", "Solon-Biet et al. 2014")

# List to hold all of the models
models<-list()
letters<-c("A", "B")

# Fit the leverage model through the isocaloic datasets
for(i in 1:2){

	# Pull out the right dataset
	data<-datasets[[i]]
	# Constrain dataset to those less than 18kJ/g
	data<-data[which(data$kJ.g < 18 & data$kJ.g > 15),]
	
	# Plot the data
	plot(data$diet.g.P, data$intake.gram, xlab="Prop. Protein by Mass", ylab="Food Intake (g)", xlim=c(0, 1), ylim=c(0, 6), pch=pcars, cex=psize)
	mtext(titles[i], line=-2, cex=1.4)
	mtext(letters[i], at=-0.05, line=1.25, cex=2, font=2)
	
	# Estimate the parameters for leverage from the data by least squares
	p<-data$diet.g.P
	intakes<-data$intake.gram
	leverage<-nls(intakes ~ P*p^L, start=list(P=2, L=-1))
	parameters<-summary(leverage)$coefficient[,1]
	
	# Check it out
	print(summary(leverage))
	
	# Fit the line - use a dashed line outside of actual data
	x<-seq(min(p), max(p), 0.001)
	y<-parameters[1] * x ^ parameters[2]
	lines(x, y, lwd=2)
	x<-seq(0, min(p), 0.001)
	y<-parameters[1] * x ^ parameters[2]
	lines(x, y, lty=2, lwd=2)
	x<-seq(max(p), 1, 0.001)
	y<-parameters[1] * x ^ parameters[2]
	lines(x, y, lty=2, lwd=2)
	
	
	# Save the model
	models[[i]]<-leverage
}

# Fit the leverage model through the full datasets
letters<-c("C", "D")
for(i in 1:2){

	# Pull out the right dataset
	data<-datasets[[i]]
	
	# Plot the data
	plot(data$diet.g.P, data$intake.gram, xlab="Prop. Protein by Mass", ylab="Food Intake (g)", xlim=c(0, 1), ylim=c(0, 6), pch=pcars, cex=psize)
	mtext(titles[i], line=-2, cex=1.4)
	mtext(letters[i], at=-0.05, line=1.25, cex=2, font=2)
	
	# Estimate the parameters for leverage from the data by least squares
	p<-data$diet.g.P
	intakes<-data$intake.gram
	leverage<-nls(intakes ~ P*p^L, start=list(P=2, L=-1))
	parameters<-summary(leverage)$coefficient[,1]
	
	# Check it out
	print(summary(leverage))
	
	# Fit the line - use a dashed line outside of actual data
	x<-seq(min(p), max(p), 0.001)
	y<-parameters[1] * x ^ parameters[2]
	lines(x, y, lwd=2)
	x<-seq(0, min(p), 0.001)
	y<-parameters[1] * x ^ parameters[2]
	lines(x, y, lty=2, lwd=2)
	x<-seq(max(p), 1, 0.001)
	y<-parameters[1] * x ^ parameters[2]
	lines(x, y, lty=2, lwd=2)
	
	# Save the model
	models[[i+2]]<-leverage

}

dev.off()

# Table to hold results
results<-as.data.frame(array(NA, c(8,7)))
names(results)<-c("Dataset", "Study", "Coef.", "Est.", "SE", "t", "p")
results$Dataset<-c(rep("15.9-18kJ/g", 4), rep("All Data", 4))
results$Study<-rep(c(rep("Hu et al.", 2), rep("Solon-Biet et al.", 2)), 2)
results$Coef.<-rep(c("P", "L"), 4)

# Order of models in table
table.order<-c(1,3,2,4)
counter<-1

# Add in the model results
for(i in 1:4){
	
	# Pull out the model and save in the table then remove it
	model<-models[[table.order[i]]]
	results[c(counter:(counter+1)),c(4:7)]<-summary(model)$coefficients
	counter<-counter+2	
	rm(model)
	
}

# Round the figures and write the results
results[,c(4:6)]<-round(results[,c(4:6)], 3)
write.table(results, file="TableS5.csv", sep=",", row.names=F, col.names=names(results))


############################################################
######################## Figure 5 ##########################
############################################################

rm(list=ls())

library(plyr)

# For this analysis we need to calculate the ratio of dry weight protein to non-protein, where non-protein is either non-energy content of the diet or non-protein energy grams, as well as energy intake

#####

# Read the animal data
data<-read.csv("Solon_Biet_et_al_Intakes.csv")
head(data)

# Read in the dietary data
diets<-read.csv("Solon_Biet_et_al_Diets.csv")

# match up the diets to the intakes and find the proportion of each macro, dry weight per gram of each animals diet
data$diet.P<-diets$diet.P[match(data$diet, diets$diet)]
data$diet.C<-diets$diet.C[match(data$diet, diets$diet)]
data$diet.F<-diets$diet.F[match(data$diet, diets$diet)]

# match up the diets to the intakes and find the proportion of each macro, dry weight per gram of each animals diet
data$diet.g.P<-diets$Dry.P[match(data$diet, diets$diet)]
data$diet.g.C<-diets$Dry.C[match(data$diet, diets$diet)]
data$diet.g.F<-diets$Dry.F[match(data$diet, diets$diet)]

# Non energy content (assumings it's cellulose - but its actually probably a bit more than that)
data$cell.g<-(1 - (data$diet.g.P + data$diet.g.C + data$diet.g.F))

# Calculate the ratio of P to cell in g, and P to fat in g
data$P.cell<-(data$diet.g.P / data$cell.g)
data$P.NP<-(data$diet.g.P / data$diet.g.F)

# match up the diets to the intakes and find the Kjs of each macro per gram of each animals diet
data$diet.P<-diets$diet.P[match(data$diet, diets$diet)]
data$diet.C<-diets$diet.C[match(data$diet, diets$diet)]
data$diet.F<-diets$diet.F[match(data$diet, diets$diet)]

# Calculate the kJs per gram
data$kJ.g<-data$diet.P + data$diet.C + data$diet.F

# Calculate energy intake
data$intake.E<-data$intake.gram * data$kJ.g

# Plotting  

pdf("Figure5.pdf", width=10, height=10)

# Aesthetics
par(mfrow=c(2,2), mar=c(5,5,5,1), cex.axis=1.2, cex.lab=1.5)
point.type<-16
point.size<-0.6

# Ratios to plot
ratios<-c("P.cell", "P.NP")
xlabs<-c("Protein (g) / Cellulose (g)", "Protein (g) / Fat (g)")

# list to hold models
models<-list()
letters<-c("A", "B")

for(i in 1:2){
	
	# Use the ith ratio as the predictor
	ith.ratio<-data[,ratios[i]]
	
	# Plot the data
	plot(ith.ratio, data$intake.gram, pch=point.type, cex=point.size, xlab=xlabs[i], ylab="Food Intake (g)", col=1, xlim=c(0, 22))
	
	# Fit the linear model, with the ratio log transformed
	model<-lm(data$intake.gram ~ log(ith.ratio))
	
	# Add the model predictions
	x<-seq(min(ith.ratio), max(ith.ratio), 0.01)
	y<-model$coef[1] + model$coef[2]*log(x)
	lines(x, y, lwd=2)
	x<-seq(max(ith.ratio), 23, 0.01)
	y<-model$coef[1] + model$coef[2]*log(x)
	lines(x, y, lwd=2, lty=2)
	
	# Save thee model
	models[[i]]<-model
	mtext(letters[i], line=1.5, at=-0.15, font=2, cex=2)
	rm(model)
	
}

letters<-c("C", "D")
for(i in 1:2){
	
	ith.ratio<-data[,ratios[i]]
	
	plot(ith.ratio, data$intake.E, pch=point.type, cex=point.size, xlab=xlabs[i], ylab="Energy Intake (kJ)", col=1, xlim=c(0, 22))
	
	# Fit the linear model, with the ratio log transformed
	model<-lm(data$intake.E ~ log(ith.ratio))
	
	# Add the model predictions
	x<-seq(min(ith.ratio), max(ith.ratio), 0.01)
	y<-model$coef[1] + model$coef[2]*log(x)
	lines(x, y, lwd=2)
	x<-seq(max(ith.ratio), 23, 0.01)
	y<-model$coef[1] + model$coef[2]*log(x)
	lines(x, y, lwd=2, lty=2)
	
	models[[i+2]]<-model
	mtext(letters[i], line=1.5, at=-0.15, font=2, cex=2)
	rm(model)
}

dev.off()

# Table to hold results
results<-as.data.frame(array(NA, c(8,7)))
names(results)<-c("Outcome", "Ratio", "Coef.", "Est.", "SE", "t", "p")
results$Outcome<-c(rep("Mass Intake", 4), rep("Energy Intake", 4))
results$Ratio<-rep(c(rep("P / cellulose", 2), rep("P / F", 2)), 2)
results$Coef.<-rep(c("Intercept", "ln Protein Ratio"), 4)

# Order of models in table
table.order<-c(1,2,3,4)
counter<-1

# Add in the model results
for(i in 1:4){
	
	# Pull out the model and save in the table then remove it
	model<-models[[table.order[i]]]
	results[c(counter:(counter+1)),c(4:7)]<-summary(model)$coefficients
	counter<-counter+2	
	rm(model)
	
}

# Round the figures
results[,c(4:6)]<-round(results[,c(4:6)], 3)
write.table(results, file="TableS6.csv", sep=",", row.names=F, col.names=names(results))


############################################################
######################## Figure S1 ##########################
############################################################

# Clear any old objects
rm(list=ls())

# Load the library for plyr
library(plyr)

############## Data from speakman paper

# Read the animal data from Hu et al.
data<-read.csv("Hu_et_al_Intakes.csv")
head(data)

# We are only interested in diets 1 through 12
data<-data[-which(data$diet > 12),]

# Read in the dietary data
diets<-read.csv("Hu_et_al_Diets_in_S1.csv")

# Convert the calories per gram to kJs
diets$kJ.g<-diets$energy.density * 4.184

# Match up the diets to the intakes and find the compositon and energy density of each animals intake
data$kJ.g<-diets$kJ.g[match(data$diet, diets$diet)]
data$diet.P.g<-diets$Dry.P[match(data$diet, diets$diet)]
data$per.P<-diets$percent.P[match(data$diet, diets$diet)]

# Calculate the values from their diet compositions
data$intake.E<-data$kJ.g * data$intake.gram
data$intake.P.g<-(data$diet.P.g/100) * data$intake.gram

# Based on Fig 2, looks like they set any NA data for protein in grams to 0, but not for Energy and overall food intake
data$intake.P.g[which(is.na(data$intake.P.g) == T)]<-0

# The data sets list
datasets<-list()
datasets[[1]]<-data

##########################################

# now recalculate these values from our derived diet data

# Read in the dietary data
diets<-read.csv("Hu_et_al_Diets.csv")

# Convert the calories per gram to kJs
diets$kJ.g<-diets$energy.density * 4.184

# Match up the diets to the intakes and find the compositon and energy density of each animals intake
data$kJ.g<-diets$kJ.g[match(data$diet, diets$diet)]
data$diet.P.g<-diets$Dry.P[match(data$diet, diets$diet)]
data$per.P<-diets$percent.P[match(data$diet, diets$diet)]

# Calculate the values from rederived diet compositions
data$intake.E<-data$kJ.g * data$intake.gram
data$intake.P.g<-data$diet.P.g * data$intake.gram

# Here we will not set any missing intakes to 0 

# Save the dataset
datasets[[2]]<-data
rm(data)

##########################################

pdf("FigureS1.pdf", height=5, width=10)

# Define the sixe of the points
p.size<-1.2

# Specify the outcomes to plot and the aesthetics
outcomes<-c("intake.E", "intake.P.g", "intake.gram")
cars<-list()
cars[[1]]<-c(18, 15, 17)
cars[[2]]<-c(5, 0, 2)
cols<-c("blue", "red", "green")

# Set up the plot region
par(mfrow=c(1,2), mar=c(5,5,5,5))

# Make empty plot
plot(-10, -10, ylim=c(0, 100), xlab="Protein content (%)", ylab="kJ/day", main="C57BL/6:60% fat", xlim=c(0, 35), pch=16, cex=p.size, bty="n", xaxs="i", yaxs="i", axes=F)
legend(20, 103, c("Energy intake", "Protein intake", "Food intake"), pch=cars[[1]], col=cols, bty="n", cex=0.8)
axis(1)
axis(2, at=seq(0, 100, 10))
axis(4, at=seq(0, 100, 10), labels=seq(0, 5, 0.5))
mtext("g/day", side=4, line=2.5, cex=1)

# Open the loop for each dataset
for(k in 1:2){
	
	# Use the kth data 1 - original diets, 2 - re-dervied diets
	data<-datasets[[k]]
	
	# Series one diets with 60% fat
	data1<-data[which(data$diet <= 6),]
	
	# Offset the points a little
	offset<-1
	if(k==1){offset<--1}
	
	# Open the loop
	for(i in 1:3){
		
		# Find the right outcome
		data1$ith.outcome<-data1[,outcomes[i]]
			
		# Scale for the mass intakes
		if(i > 1){
			data1$ith.outcome<-(data1$ith.outcome/5) * 100
		}
		
		# Calculate the mean and SD
		means<-ddply(data1, .(diet), summarise, mean(ith.outcome, na.rm=T), sd(ith.outcome, na.rm=T), per.P=mean(per.P))
		
		# Add in the points and SDs
		arrows(means$per.P+offset, means$..1 + means$..2, means$per.P+offset, means$..1 - means$..2, code=3, angle=90, len=0.02)
		
		# If on K2 - overlay a white point first to cover the arrows
		if(k == 2){points(means$per.P+offset, means$..1, col="white", pch=cars[[1]][i], cex=p.size*1.3)}
		points(means$per.P+offset, means$..1, col=cols[i], pch=cars[[k]][i], cex=p.size)	
	}
}

# Re plot for the second series of diets

# Make empty plot
plot(-10, -10, ylim=c(0, 100), xlab="Protein content (%)", ylab="kJ/day", main="C57BL/6:20% fat", xlim=c(0, 35), pch=16, cex=p.size, bty="n", xaxs="i", yaxs="i", axes=F)
axis(1)
axis(2, at=seq(0, 100, 10))
axis(4, at=seq(0, 100, 10), labels=seq(0, 5, 0.5))
mtext("g/day", side=4, line=2.5, cex=1)

# Open the loop for the two datasets
for(k in 1:2){
	
	# Use the kth data 1 - original diets, 2 - re-dervied diets
	data<-datasets[[k]]
	
	# Series two diets with 20% fat
	data1<-data[which(data$diet > 6),]
	
	# Offset the points a little
	offset<-1
	if(k==1){offset<--1}
	
	# Open the loop
	for(i in 1:3){
		
		# Find the right outcome
		data1$ith.outcome<-data1[,outcomes[i]]
			
		# Scale for the mass intakes
		if(i > 1){
			data1$ith.outcome<-(data1$ith.outcome/5) * 100
		}
		
		# Calculate the mean and SD
		means<-ddply(data1, .(diet), summarise, mean(ith.outcome, na.rm=T), sd(ith.outcome, na.rm=T), per.P=mean(per.P))
		
		# Add in the points and SDs
		arrows(means$per.P+offset, means$..1 + means$..2, means$per.P+offset, means$..1 - means$..2, code=3, angle=90, len=0.02)
		
		# If on K2 - overlay a white point first to cover the arrows
		if(k == 2){points(means$per.P+offset, means$..1, col="white", pch=cars[[1]][i], cex=p.size*1.2)}
		points(means$per.P+offset, means$..1, col=cols[i], pch=cars[[k]][i], cex=p.size)
		
	}

}

dev.off()
