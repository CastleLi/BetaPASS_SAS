*=========================================  Macro part  =================================================*;
*========================================================================================================*;

**SAS program for calculating power and sample size for the beta distribution 6-10-15 ;
options nocenter formdlim=' ' pagesize=100 linesize=128; run;
 OPTIONS FORMCHAR="|----|+|---+=|-/\<>*"; run;
 ods listing;
 ods html close; run;

options nomacrogen;
options nonotes;

** return power with given sample size and link type;
%macro betapwr2(sampsize);
data powerest;powerhat=.;run;

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


%if &link_type ^= wilcoxon %then %do;
*=================== calculate power of GLIMMIX method ==================*;
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

data powerest;
  set outm2(keep = power);
  rename power = powerhat;
run;




%end;
%else %do;
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
output out=outw2 mean=power;
title2 "Power with mu1=&mu1, &sampsize sampsize, seed=&seed"; run;

data powerest;
  set outw2(keep = power);
  rename power = powerhat;
run;
%end;

%mend; run;



* return sample size with given mu1, target power and link type;
%macro sample_size_mid(mu1,power_target,delta,seed,link_type);

ods select none;
proc power;
      twosamplemeans
         nfractional
         meandiff     = %sysevalf(&mu1 - &mu0)
         stddev       = &sd0
         power        = &power_target
         ntotal       = .;
		 ODS Output Output = startpw;
		 run;
ods select all;
data _null_;
  set startpw;
  call symput("ss_start",NTotal);
run;


%betapwr2(&ss_start);
data _null_;
set powerest;
call symput("powerhat",powerhat);
run;

%if %sysevalf(&powerhat) < &power_target %then %do;
	%let ss_end = %sysevalf(2*&ss_start,floor);
	%betapwr2(&ss_end);
	data _null_;
 	set powerest;
    call symput("powerhat",powerhat);
	run;
	%do %while(&powerhat < &power_target);
		%let ss_start = %sysevalf(&ss_end);
		%let ss_end = %sysevalf(2*&ss_start,floor);
		%betapwr2(&ss_end);
		data _null_;
		set powerest;
		call symput("powerhat",powerhat);
		run;
	%end;
	%let ss_lower = &ss_start;
	%let ss_upper = &ss_end;
%end;

%else %do;
	%let ss_end = %sysevalf(0.5*&ss_start,floor);
	%betapwr2(&ss_end);
	data _null_;
    set powerest;
    call symput("powerhat",powerhat);
	run;
	%do %while((&powerhat > &power_target) and (&ss_end >5));
		%let ss_start = %sysevalf(&ss_end);
		%let ss_end = %sysevalf(0.5*&ss_start,floor);
		%betapwr2(&ss_end);
		data _null_;
 		set powerest;
 		call symput("powerhat",powerhat);
		run;
	%end;
	%let ss_lower = &ss_end;
	%let ss_upper = &ss_start;
%end;

%put &ss_lower;
%put &ss_upper;

%let reachflag =0;
%do %while(&reachflag<1);
	%let ss_diff = %sysevalf(&ss_upper - &ss_lower,int);
	%if &ss_diff <= &delta  %then %let reachflag = 1;
	%else %do;
		%let ss_mid = %sysevalf((&ss_upper + &ss_lower)*0.5,floor);
		%betapwr2(&ss_mid);
		data _null_;
 		set powerest;
 		call symput("powerhat",powerhat);
		run;
		%if &powerhat < &power_target %then %do;
			%let ss_lower = &ss_mid;
		%end;
		%else %do;
			%let ss_upper = &ss_mid;
		%end;
	%end;

%put &ss_lower;
%put &ss_upper;

%end;

%let ss_min = &ss_upper;
%betapwr2(&ss_min);
data _null_;
set powerest;
call symput("powerhat",powerhat);
run;


data ssest;
mu1 = &mu1;
minss = &ss_min;
minpower = &powerhat;
run;
%let minss = minss;
%let minpower = minpower;
data ssest;
set ssest;
rename &minss = minimum_samplesize_of_&link_type;
rename &minpower = minimum_power_of_link_&link_type;
run;

%mend;run;



%macro doit_ss1;
data oness; 
mu1 = &mu1;
target_power = &power_target;
run;

%if %scan(&_link_type, 1) = all %then %do;
%let _link_type = logit probit log cloglog loglog;
%end;
%let nlinktype=%sysfunc(countw(&_link_type));
%do l=1 %to &nlinktype;
%let link_type = %scan(&_link_type, &l);

%sample_size_mid(&mu1,&power_target,&delta,&seed,&link_type);
data oness;
merge oness ssest;
by mu1;
run;
%end;

data allss;
set allss oness;
run;


%mend;run;


%macro doit_ss2;
%if &mu1_by = 0 %then %do;
%let loops = 0;
%end;
%else %do;
%let loops=%sysevalf((&mu1_end-&mu1_start)/&mu1_by); 
%end;
%if &power_by = 0 %then %do;
%let loops2 = 0;
%end;
%else %do; 
%let loops2 = %sysevalf((&power_end-&power_start)/&power_by); 
%end;

%do j=0 %to &loops2;
  %let power_target = %sysevalf(&power_start+&j*&power_by);  
  %let mu1=%sysevalf(&mu1_start); 
  %do i=0 %to &loops;
    %let seed=%eval (&seed +1);
    %doit_ss1
    %let mu1=%sysevalf(&mu1+&mu1_by);  
	%put target power = &power_target;
    %put mu1= &mu1; 
    %put seed=&seed;
  %end;
%end;

%mend; run;

%macro samplesize(mu0,sd0,mu1_start,mu1_end,mu1_by,power_start,power_end,power_by,delta,trials,seed,_link_type,equal_precision,sd1);
data allss;ini=.;
%doit_ss2;

data allss;
set allss(drop = ini);
run;
data &output_matrix_ss;
set allss;
where mu1^=.;
run;

proc print data=&output_matrix_ss; 
title2 "based on &trials trials";
run;
%mend;run;


%macro plot_samplesize(_link_type);

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
	set &output_matrix_ss;
	rename minimum_samplesize_of_&link_type = samplesize;
	rename minimum_power_of_link_&link_type = power;
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

proc sort data=allmt; 
by mu1;
run;

proc sgpanel data=allmt; 
panelby mu1;
series x=power y=samplesize/ markers group=method;
title "Minimum power and sample size, by mu1";
run;

%mend;run;

%let output_matrix_ss = Test1;
%let _link_type = log logit;
%samplesize(0.56,0.255,.70,.75,.05,0.7,0.9,0.1,1,40,610201501,&_link_type,equal_precision = TRUE,sd1 = 0.01); run;
%samplesize(0.56,0.255,.70,.75,.05,0.7,0.9,0.1,1,40,610201501,&_link_type,equal_precision = FALSE,sd1 = 0.01); run;
proc print data=Test1;run;
%plot_samplesize(&_link_type);

%let output_matrix_ss = Test2;
%let _link_type = all;
%samplesize(0.56,0.255,.65,.75,.05,0.7,0.9,0.1,1,200,610201501,&_link_type,equal_precision = TRUE,sd1 = 0.01); run;
%samplesize(0.56,0.255,.65,.75,.05,0.7,0.9,0.1,1,200,610201501,&_link_type,equal_precision = FALSE,sd1 = 0.01); run;
proc print data=Test2;run;
%plot_samplesize(&_link_type);

*=========================================  Example part  ===============================================*;
*========================================================================================================*;

%let output_matrix_ss = Example1;
%let _link_type = logit wilcoxon;
%samplesize(0.0174,0.0211,.0131,.0131,0,0.8,0.8,0,1,1000,1,&_link_type,TRUE); run;
%samplesize(0.0174,0.0211,.0131,.0131,0,0.8,0.8,0,1,1000,1,&_link_type,FALSE,0.01); run;
