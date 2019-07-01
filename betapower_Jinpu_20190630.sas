*=========================================  Macro part  =================================================*;
*========================================================================================================*;

**SAS program for calculating power and sample size for the beta distribution 6-10-15 ;
options nocenter formdlim=' ' pagesize=100 linesize=128; run;
 OPTIONS FORMCHAR="|----|+|---+=|-/\<>*"; run;
 ods listing;
 ods html close; run;

**Assume the response variable Y is a proportion or otherwise constrained between 0 and 1 which can be modeled with a beta dbn.
** This program will help find the power for a given sample size when testing the null hypothesis that the means for the
** control and treatment groups are equal against a two-sided alternative.
** The user must supply the mean and standard deviation for the control group (mu0 and sd0) as well as the mean for the
** treatment group under the alternative, namely mu1.
** See SAS documentation under DIST=BETA as a model option for GLIMMIX.  If Y~beta(a,b), then mu=a/(a+b) and the variance of Y
** can be expressed in terms of mu using the parameter phi as Var(Y)=mu*(1-mu)/(1+phi).  The value of phi is found from mu0
** and sd0.  That value is then used to find the variance under the alternative.  Given mu and phi the parameters a and b
** can be found from a=mu*phi and b=(1-mu)*phi.  The values of a and b are then used to generate random beta variables for
** the simulation.
 *****************************;
** The macro allows you to control the number of trials in the simulation, the sample sizes used, and the alternative means.
** You can fix the alternative and vary sample size to match a desired power;
** You can fix the sample size and vary the alternative to see which will match a desired power;
** You can vary both;
  **Start with a small number of trials (say 100) to determine the rough range of sample sizes or alternatives;
  ** Use a larger number of trials (say 1000) to get better estimates;
 ** The output includes SGPLOTS of Power curves;
 **
 ** The last macro (BETAPOWER) is the one to which you supply the parameters.
 ** The betapower macro invokes the DOIT macro which in turn calls the BETAPWR macro;
*********;
options nomacrogen;
options nonotes;

%macro betapwr(mu1,sampsize,seed);

*=================== generate data ==================*;
data sim;  
*** set seed;
call streaminit(&seed);
*** input mu0 and mu1;
mu0=&mu0; 
mu1=&mu1;
*** set parameters;
%if &equal_precision = TRUE %then %do;
phi=((mu0*(1-mu0))/(&sd0*&sd0))-1;
a0=mu0*phi; 
b0=(1-mu0)*phi;
a1=mu1*phi; 
b1=(1-mu1)*phi;
%end;
%else %do;
phi=((mu0*(1-mu0))/(&sd0*&sd0))-1;
phi1=((mu1*(1-mu1))/(&sd1*&sd1))-1;
a0=mu0*phi; 
b0=(1-mu0)*phi;
a1=mu1*phi1; 
b1=(1-mu1)*phi1;
%end;
*** set dependent variable;
do trial=1 to &trials;
  do sample = 1 to &sampsize;
    tmt=0;
    y=rand('beta',a0,b0); * y=max(y,0.00001);   
	output;  
    tmt=1; 
	y=rand('beta',a1,b1); * y=max(y,0.00001);   
	output;  
  end; 
end; 
run;


*=================== calculate power of GLIMMIX method ==================*;
data outm1;ini=1;


%if %scan(&_link_type, 1) = all %then %do;
%let _link_type = logit probit log cloglog loglog;
%end;

%let nlinktype=%sysfunc(countw(&_link_type));
%do l=1 %to &nlinktype;
%let link_type = %scan(&_link_type, &l);

%if &equal_precision = TRUE %then %do;
proc glimmix data=sim; 
by trial;  * where trial<=2;
ods output Tests3=outtest;
ods listing select none;
 class  tmt;
 model y= tmt   / s cl ddfm=kr  dist=beta link=&link_type;
 * lsmeans tmt/ilink diff cl;
 * output out=resid1 pearson=pearson  student=student;
  title "Generalized linear model: beta, &link_type link";
 run;
 ods listing select all;

data outtest; 
set outtest;
if probf^=. then  p05=(probf<0.05);

proc sort data=outtest; 
by effect;
proc means data=outtest noprint; 
by effect; 
var  p05;
output out=outm2 mean=power;
run;
%end;

%else %do;
proc nlmixed data =sim;
by trial;  * where trial<=2;
ods output ParameterEstimates=outtest;
ods listing select none;
parameters b0 = -1
b1 = -0.1
d0 = -1;
xb = b0+b1*tmt;

%if &link_type = logit %then %do;
    mu = exp(xb)/(1 + exp(xb)); 
	%put "link type = " &link_type;
%end; * logit link;
%if &link_type = probit %then %do;
      mu = cdf('NORMAL',xb); * probit link;
	  %put "link type = " &link_type;
	%end;
%if &link_type = log %then %do;
        mu = exp(xb); * log link;
		%put "link type = " &link_type;
      %end;
%if &link_type = cloglog %then %do;
          mu = 1-exp(-exp(xb)); * complementary log-log link;
		  %put "link type = " &link_type;
		%end;
%if &link_type = loglog %then %do;
            mu = exp(-exp(-xb)); * log-log link;
			%put "link type = " &link_type;
		  %end;


phi = exp(d0);
 w = mu*phi;
 t = phi - mu*phi;
 ll = lgamma(w+t) - lgamma(w) - lgamma(t) +
 ((w-1)*log(y)) + ((t-1)*log(1 - y));
 model y ~ general(ll);
run; 
 ods listing select all;
data outtest; 
set outtest;
where Parameter = "b1"; 
run;
data outtest;
set outtest;
if Probt~=. then  p05=(Probt<0.05);
proc means data=outtest noprint; 
var  p05;
output out=outm2 mean=power;
run;
%end;

data outm2;
set outm2(keep=power);
rename power = power_mid;
run;
%let power_mid = power_mid;
data outm2;
set outm2;
ini =1;
rename &power_mid = power_of_link_&link_type;

data outm1;
merge outm1 outm2;
by ini;
run;
%end;

data outm1; 
set outm1(drop = ini); 
samplesize=&sampsize; 
mu1=&mu1;
data outm1;
set outm1;
data allm1; 
set allm1 
outm1;

*=================== calculate power of Wilcoxon Rank sum test ==================*;
proc npar1way data=sim wilcoxon; by trial; * where trial<=2;
var y; 
class tmt;
ods listing select none;
ods output "Two-Sample Test"=outwilc;
run;
ods listing select all;

data outwilc; 
set outwilc; 
if name1="PT2_WIL";
p05=(nvalue1<0.05);

proc sort data=outwilc; 
by name1;
proc means data=outwilc noprint; 
by name1; var  p05;
output out=outw1 mean=power;
title2 "Power with mu1=&mu1, &sampsize sampsize, seed=&seed"; run;
data outw1; 
set outw1(keep =power);  
samplesize=&sampsize; 
mu1=&mu1;
rename power = power_of_Wilcoxon_test;
data allw1; 
set allw1 outw1;

%mend; run;



%macro doit;
%if &mu1_by = 0 %then %do;
%let loops = 0;
%end;
%else %do;
%let loops=%sysevalf((&mu1_end-&mu1_start)/&mu1_by); 
%end;

%do sampsize=&ss_start %to &ss_end %by &ss_by;
  %let mu1=%sysevalf(&mu1_start); 
  %put first mu1= &mu1; 
  %put mu1_end=&mu1_end; 
  %put mu1_by= &mu1_by;
  %do i=0 %to &loops;
    %let seed=%eval (&seed +1);
    %betapwr(&mu1,&sampsize,&seed); 
    %let mu1=%sysevalf(&mu1+&mu1_by);  
    %put second mu1= &mu1; 
    %put seed=&seed;
  %end;
%end;

%mend; run;








%macro betapower(mu0,sd0,mu1_start,mu1_end,mu1_by,ss_start,ss_end,ss_by,trials,seed,_link_type,equal_precision,sd1);
data allm1; power=.;
data allw1; power=.; 

%doit;

data allmw(drop = power);
merge allm1 allw1(in = in2);
by samplesize mu1;
run;

data &output_matrix;
set allmw;
where samplesize^=.;
run;

proc print data=&output_matrix; 
title2 "based on &trials trials";
run;

%mend; run;


%macro plot_betapower(Plotby,_link_type);



data allmt;
power = .;
run;

%if %scan(&_link_type, 1) = all %then %do;
%let _link_type = logit probit log cloglog loglog;
%end;

%let nlinktype=%sysfunc(countw(&_link_type));
%do l=1 %to &nlinktype;
%let link_type = %scan(&_link_type, &l);
	data onemt;
	set &output_matrix;
	rename power_of_link_&link_type = power;
	run;
	data onemt;
	set onemt(keep = power samplesize mu1);
	method= "&link_type";
	run;

	data allmt;
	length method $17;
	set allmt onemt;
	where power^=.;
	run;
%end;

data wlmt;
set  &output_matrix(keep = power_of_Wilcoxon_test samplesize mu1);
rename power_of_Wilcoxon_test=power;
method="WILCOXON";
run;

data allmw; 
set allmt wlmt;
run;

proc sort data=allmw; 
by mu1;
run;

%if &plotby = linktype %then %do;
	proc sgpanel data=allmw;
	panelby method;
	series x=mu1 y=power/ markers group=samplesize;
	refline 0.80;
	title "Power for GLIMMIX and Wilcoxon test, by link type";
	run;
%end;
%else %do;
	%if &plotby = samplesize %then %do;
		proc sgpanel data=allmw; 
		panelby samplesize;
		series x=mu1 y=power/ markers group=method;
		refline 0.80;
		title "Power for GLIMMIX and Wilcoxon test, by sample size";
		run;
	%end;

	%else %do;
		proc sgpanel data=allmw; 
		panelby mu1;
		series x=samplesize y=power/ markers group=method;
		refline 0.80;
		title "Power for GLIMMIX and Wilcoxon test, by mu1";
		run;
	%end;
%end;

%mend; run;







*=========================================  Test part  ==================================================*;
*========================================================================================================*;


*betapower(mu0,sd0,mu1_start,mu1_end,mu1_by,ss_start,ss_end,ss_by,trials,seed);
%let output_matrix = Test1;
%let _link_type = logit log;
%betapower(0.56,0.255,.70,.75,.05,30,50, 20,40,1,&_link_type,TRUE); run;
%betapower(0.56,0.255,.70,.75,.05,30,50, 20,40,1,&_link_type,FALSE,0.2); run;
proc print data=Test1;run;
%plot_betapower(linktype,&_link_type);run;
%plot_betapower(samplesize,&_link_type);run;
%plot_betapower(mu1,&_link_type);run;

%let output_matrix = Test2;
%let _link_type = logit log;
%betapower(0.56,0.255,.60,.75,.05,30,50, 5,1000,1,&_link_type,TRUE,0.01); run;
%betapower(0.56,0.255,.60,.75,.05,30,50, 5,1000,1,&_link_type,FALSE,0.01); run;
proc print data=Test2;run;
%plot_betapower(linktype,&_link_type);run;
%plot_betapower(samplesize,&_link_type);run;
%plot_betapower(mu1,&_link_type);run;

%let output_matrix = Test3;
%let _link_type = all;
%betapower(0.56,0.255,.560,.76,.05,30,45, 5,200,610201511,&_link_type,TRUE,0.2); run;
%betapower(0.56,0.255,.560,.76,.05,30,45, 5,200,610201511,&_link_type,FALSE,0.2); run;
%plot_betapower(linktype,&_link_type);run;
%plot_betapower(samplesize,&_link_type);run;
%plot_betapower(mu1,&_link_type);run;


%let output_matrix = Test4;
%let _link_type = logit;
%betapower(0.56,0.255,.560,.56,.05,30,45, 5,100,617201521,&_link_type,TRUE,0.01); run;
%betapower(0.56,0.255,.560,.56,.05,30,45, 5,100,617201521,&_link_type,FALSE,0.01); run;
%plot_betapower(linktype,&_link_type);run;
%plot_betapower(samplesize,&_link_type);run;
%plot_betapower(mu1,&_link_type);run;





*=========================================  Example part  ===============================================*;
*========================================================================================================*;

%let output_matrix = Example1;
%let _link_type = logit;
%betapower(0.0174,0.0211,0.0120,0.0140, 0.0010,100,200, 25,500,1,&_link_type,TRUE,0.01); run;
%plot_betapower(linktype,&_link_type);run;
%betapower(0.0174,0.0211,0.0120,0.0140, 0.0010,100,200, 25,500,1,&_link_type,FALSE,0.01); run;
%plot_betapower(linktype,&_link_type);run;
