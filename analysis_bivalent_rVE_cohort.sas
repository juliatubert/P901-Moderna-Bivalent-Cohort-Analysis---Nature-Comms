/*********************************************************************************************************************************************************
Program: analysis_bivalent_rVE_cohort.sas

Purpose: 
    Moderna bivalent cohort analysis:
	For rVE cohort: Bivalent-vaccinated vs. individuals with 2+ mRNA doses only
    Describe the study population baseline characteristics  
    and conduct the analysis for study outcomes.

Note: Table and figure #s may differ in final publication 
*********************************************************************************************************************************************************/

*** Set up infile passway;
libname in "~/bivalent_formal" access=readonly;
%let infile = in.covariates_cohort_rVE_CR;

*** Set up output passway;
%let outpath = bivalent_analysis_formal;
libname out  "~/analysis/&outpath";
%let update_date=%sysfunc(today(), date9.); 
%put &update_date;

*** Follow up;
%let fu_end='31Jan2023'd;

*** Study arms;
%let exp=Bivalent;
%let unexp=TwoDose;

*** Call macros;
%include "~/analysis/&outpath/macros_bivalent.sas";


**********************************
** Table 1. Baseline characteristics;

* Prepare data;
    proc univariate data=&infile plot; 
    var Frailty_index;
    ods output Quantiles=Qs;
    run;
    proc sql; 
    select estimate into: Q1 from Qs where quantile='25% Q1';
    select estimate into: Q2 from Qs where quantile='50% Median';
    select estimate into: Q3 from Qs where quantile='75% Q3';
    quit;
    data data_baseline; set &infile;
    * Create frailty category based on top qualtile;
        if Frailty_index >=&Q3 then Frailty=4;
        else if Frailty_index >=&Q2 then Frailty=3;
        else if Frailty_index >=&Q1 then Frailty=2;
        else Frailty=1;
    * Adjust variable format and create new variables for stddiff macro;
        ARM_n=(ARM="&exp"); 
        sex_n=(sex='Male');
    run; 
	* Remove variable labels to simplify baseline table generation;
	proc datasets lib=work;
		modify data_baseline;
		attrib _all_ label='';
	run;
	proc freq data=data_baseline;
	table Frailty ARM*ARM_n*cohort sex*sex_n/list missing; 
	run; 

* Run baseline table macro;
	%baselinetable(data_baseline, table1, ARM_n);
	title 'Table 1';
	Proc print data=table1; run;
	title;

* Check covariates with ASD>0.1 for model adjustment;
	proc sql;
		create table significant as
		select a.VarName, a.Stddiff_4dec as ABS_4dec, b.ABS as ABS_table1  
		from table1_std a
		left join table1 b on a.VarName=b.var
		where (a.Stddiff_4dec>0.1 or .<a.Stddiff_4dec<-0.1) and b.ABS^='N/A'
		order by VarName;
	quit; 
	*VA_grp days_since_hist_covid_c days_since_last_vac days_since_last_vac_c medical_center_and_area ndose_c preventive_care;
	*Look out for continuous vs categorical vars and vars that are only applicable to subgroups;
		 
* Add footnote for concomitant vaccinations;
	proc freq data=data_baseline;
		table concomitant_influenza concomitant_pneumococcal concomitant_Tdap concomitant_zoster concomitant_other_vac;
		where concomitant_vac;
	run;
	proc sql noprint; 
	   select count(*) into: concomvac from data_baseline where concomitant_vac;
	   select put(count(distinct mrn)/&concomvac*100, 4.1) into: inf  trimmed from data_baseline where concomitant_influenza;
	   select put(count(distinct mrn)/&concomvac*100, 4.1) into: zost trimmed from data_baseline where concomitant_zoster;
	   select put(count(distinct mrn)/&concomvac*100, 4.1) into: pne  trimmed from data_baseline where concomitant_pneumococcal;
	   select put(count(distinct mrn)/&concomvac*100, 4.1) into: tdap trimmed from data_baseline where concomitant_Tdap;
	   select put(count(distinct mrn)/&concomvac*100, 4.1) into: oth  trimmed from data_baseline where concomitant_other_vac;
	quit;
	%put "Among subjects with concomitant vaccines received with the Moderna bivalent vaccine: influenza vaccine (&inf%), shingles vaccine (&zost%), pneumococcal vaccine (&pne%), Tdap (&tdap%), and other vaccine (&oth%)";


**************** Analysis of primary objectives ******************

** Prepare data;
    *The person-years will be the time from 14 days after the index date to the date of the event of interest, 
      death, disenrollment (allowing for a 31-day gap), end of follow-up of analysis, or 
      receipt of a dose of any COVID-19 vaccine, whichever comes first.;

* Prepare data for analysis;
    data data_analysis; set data_baseline;
    * Calculate person year for Cox models;
		py_infec_all=(min(covid_infec_all_date, dod, memb_stop-1, &FU_end, post_vacdt)-(index_date+14)+1)/365.25; 
		py_infec_edu=(min(covid_infec_edu_date, dod, memb_stop-1, &FU_end, post_vacdt)-(index_date+14)+1)/365.25;
        py_hosp=(min(covid_hosp_adate, dod, memb_stop-1, &FU_end, post_vacdt)-(index_date+14)+1)/365.25;
        py_dead=(min(dod, memb_stop-1, &FU_end, post_vacdt)-(index_date+14)+1)/365.25;
    * Calculate person months for KM-plots;
		py_infec_all_mons=py_infec_all*12;
		py_infec_edu_mons=py_infec_edu*12;
        py_hosp_mons=py_hosp*12;
        py_dead_mons=py_dead*12;
    * Calculate log person year for incidence rate calculation in macro;
		lpy_infec_all=log(py_infec_all);
		lpy_infec_edu=log(py_infec_edu);
        lpy_hosp=log(py_hosp);
        lpy_dead=log(py_dead);
    * Specify competing risk;
        covid_infec_all_comp=covid_infec_all; if (dod-(index_date+14)+1)/365.25=py_infec_all and covid_infec_all=0 then covid_infec_all_comp=2;
        covid_infec_edu_comp=covid_infec_edu; if (dod-(index_date+14)+1)/365.25=py_infec_edu and covid_infec_edu=0 then covid_infec_edu_comp=2;
        covid_hosp_comp=covid_hosp; if (dod-(index_date+14)+1)/365.25=py_hosp and covid_hosp=0 then covid_hosp_comp=2;
        covid_hosp_death_comp=covid_hosp_death; if (dod-(index_date+14)+1)/365.25=py_dead and covid_hosp_death=0 then covid_hosp_death_comp=2;
	* Collapse vars for convergence issues;
		race_eth_comb4 = race_eth; if race_eth in ('Non-Hispanic Black','Other/Unknown') then race_eth_comb4='Other/Unknown';
		age_grp_comb4 = age_grp; if age_grp in ('6-17 yo','18-44 yo') then age_grp_comb4='6-44 yo';
		age_grp_comb3 = age_grp; if age_grp in ('6-17 yo','18-44 yo','45-64 yo') then age_grp_comb3='<65 yo';
		VA_grp_comb3 = VA_grp; if VA_grp in ('0','1-4') then VA_grp_comb3='0-4';
		ndose_c_comb2=ndose_c; if ndose>=3 then ndose_c_comb2='3+';
    run;


***************************************************************;
** Table 3. rVE against COVID-19 hospitalization overall and by subgroups;

* Hospitalization rVE overall (Collapsed age_grp and dropped med center for all hosp models);
	%let varlist=age_grp_comb4 sex race_eth index_month va_grp preventive_care ndose_c days_since_last_vac_c antiviral hist_covid_diag; 
    %cohort_analysis(data_analysis, covid_hosp, &exp., &unexp., py_hosp, lpy_hosp, py_hosp_mons, hosp_rVE, survivalplot='Y'); *change survivalplot='Y';
* Hospitalization rve by subgroups;
	* Age (No cases in 6-17yo, will get warning. Adjust for continuous age. Collapsed VA_grp);
	data data_analysis_lt75 data_analysis_75up; 
	set data_analysis; 
	if .<age<75 then output data_analysis_lt75;
	if age>=75 then output data_analysis_75up;
	run;
	%let varlist=sex race_eth index_month VA_grp_comb3 preventive_care ndose_c days_since_last_vac_c antiviral hist_covid_diag; 
	%let age_term=age;
    %cohort_analysis(data_analysis_lt75, covid_hosp, &exp., &unexp., py_hosp, lpy_hosp, py_hosp_mons, age_grp_lt75, subgroup='Y', by=age_grp);
	%let age_term=age age*age;
    %cohort_analysis(data_analysis_75up, covid_hosp, &exp., &unexp., py_hosp, lpy_hosp, py_hosp_mons, age_grp_75up, subgroup='Y', by=age_grp);
	* Sex;
	%let varlist=age_grp_comb4 sex race_eth index_month va_grp preventive_care ndose_c days_since_last_vac_c antiviral hist_covid_diag;  
    %cohort_analysis(data_analysis, covid_hosp, &exp., &unexp., py_hosp, lpy_hosp, py_hosp_mons, sex, subgroup='Y', by=sex);
	* Race/ethnicity (Other/unknown not included in results);
	%let varlist=age_grp_comb4 sex race_eth index_month va_grp preventive_care ndose_c days_since_last_vac_c antiviral hist_covid_diag; 
    %cohort_analysis(data_analysis, covid_hosp, &exp., &unexp., py_hosp, lpy_hosp, py_hosp_mons, race_eth, subgroup='Y', by=race_eth);
	* IC status (Adjusted for IC subconditions in 'Yes' group);
	data data_analysis_IC0; set data_analysis; where IC=0; run;
	%let varlist=age_grp_comb4 sex race_eth_comb4 index_month va_grp preventive_care ndose_c days_since_last_vac_c antiviral hist_covid_diag; 
    %cohort_analysis(data_analysis_IC0, covid_hosp, &exp., &unexp., py_hosp, lpy_hosp, py_hosp_mons, IC0);
	data data_analysis_IC1; set data_analysis; where IC=1; run;
	%let varlist=age_grp_comb4 sex race_eth_comb4 index_month va_grp preventive_care ndose_c days_since_last_vac_c antiviral hist_covid_diag IC_DX IC_HIV IC_MEDS IC_TX; 
    %cohort_analysis(data_analysis_IC1, covid_hosp, &exp., &unexp., py_hosp, lpy_hosp, py_hosp_mons, IC1);
	* History of COVID (Adjusted for days_since_hist_covid_c in 'Yes' group);
	data data_analysis_hist0; set data_analysis; where hist_covid_diag=0; run;
	%let varlist=age_grp_comb4 sex race_eth index_month va_grp preventive_care ndose_c days_since_last_vac_c antiviral hist_covid_diag; 
	%cohort_analysis(data_analysis_hist0, covid_hosp, &exp., &unexp., py_hosp, lpy_hosp, py_hosp_mons, HIST0);
	data data_analysis_hist1; set data_analysis; where hist_covid_diag=1; run;
	%let varlist=age_grp_comb4 sex race_eth index_month va_grp preventive_care ndose_c days_since_last_vac_c antiviral days_since_hist_covid_c; 
	%cohort_analysis(data_analysis_hist1, covid_hosp, &exp., &unexp., py_hosp, lpy_hosp, py_hosp_mons, HIST1);
	* Number of monovalent doses prior to index (Dropped days_since_last_vac_c);
	%let varlist=age_grp_comb4 sex race_eth index_month va_grp preventive_care ndose_c /*days_since_last_vac_c*/ antiviral hist_covid_diag; 
    %cohort_analysis(data_analysis, covid_hosp, &exp., &unexp., py_hosp, lpy_hosp, py_hosp_mons, ndose_c, subgroup='Y', by=ndose_c);
	* Time between latest monovalent dose and index (Dropped ndose_c);
	%let varlist=age_grp_comb4 sex race_eth index_month va_grp preventive_care /*ndose_c*/ days_since_last_vac_c antiviral hist_covid_diag; 
    %cohort_analysis(data_analysis, covid_hosp, &exp., &unexp., py_hosp, lpy_hosp, py_hosp_mons, days_since_last_vac_c, subgroup='Y', by=days_since_last_vac_c);

* Hospitalization by time since bivalent dose;
	* Max follow-up = 1/31/23-(8/31/22+14) = 139 days (~4.6 mo);
	proc sql; select max(py_hosp_mons) as max_py_hosp_mons from data_analysis; quit;

	*4 intervals x 1mo: 0-1mo, 1-2mo, 2-3mo, 3+mo;
    data data_analysis_bytime; set data_analysis;
    array py_hospn{4} py_hosp1-py_hosp4;
    array lpy_hospn{4} lpy_hosp1-lpy_hosp4;
    array covid_hospn{4} covid_hosp1-covid_hosp4;
    array covid_hosp_compn{4} covid_hosp1_comp covid_hosp2_comp covid_hosp3_comp covid_hosp4_comp;
    do i=1 to 4; 
    * Update outcome status;
        covid_hospn[i]=((i-1)<=py_hosp*12<i and covid_hosp); 
        if (i-1)<=py_hosp*12<i then covid_hosp_compn[i]=covid_hosp_comp; else covid_hosp_compn[i]=0; 
    * Allow last interval to go through end of f/u;
        covid_hosp4=((4-1)<=py_hosp*12 and covid_hosp); 
        if (4-1)<=py_hosp*12 then covid_hosp_comp4=covid_hosp_comp; else covid_hosp_comp4=0; 
    * Calculate person year for Cox models;
        py_hospn[i]=max(0,min(1/12,py_hosp-(i-1)*1/12));
    * Calculate log person year for incidence rate calculation in macro;
        lpy_hospn[i]=log(py_hospn[i]);
    end;
    run;
	proc freq data=data_analysis_bytime;
		table covid_hosp*covid_hosp1*covid_hosp2*covid_hosp3*covid_hosp4/list missing;
	run; *All events should fall under one of the time categories;

	%let varlist=age_grp_comb4 sex race_eth index_month va_grp preventive_care ndose_c days_since_last_vac_c antiviral hist_covid_diag; 
    %cohort_analysis(data_analysis_bytime, covid_hosp1, &exp., &unexp., py_hosp1, lpy_hosp1, NA, time_grp1);
    %cohort_analysis(data_analysis_bytime, covid_hosp2, &exp., &unexp., py_hosp2, lpy_hosp2, NA, time_grp2); 
    %cohort_analysis(data_analysis_bytime, covid_hosp3, &exp., &unexp., py_hosp3, lpy_hosp3, NA, time_grp3); 
	* Collapsed race/eth and VA_grp; 
	%let varlist=age_grp_comb4 sex race_eth_comb4 index_month VA_grp_comb3 preventive_care ndose_c days_since_last_vac_c antiviral hist_covid_diag; 
    %cohort_analysis(data_analysis_bytime, covid_hosp4, &exp., &unexp., py_hosp4, lpy_hosp4, NA, time_grp4); 

* Set results together;
    data blank; by=' '; run;
    data table3; 
	length by $40 description $40;
    set rslt_hosp_rVE blank
		rslt_age_grp_lt75 rslt_age_grp_75up blank 
		rslt_sex blank
		rslt_race_eth blank
		rslt_IC1 rslt_IC0 blank
		rslt_hist1 rslt_hist0 blank
		rslt_time_grp1 rslt_time_grp2 rslt_time_grp3 rslt_time_grp4 blank 
		rslt_ndose_c blank
		rslt_days_since_last_vac_c 
		;
    order=_N_;
    if by='Hispanic' then order=_N_+2;
    if by='Non-Hispanic Asian' then order=_N_+2;
    if by='Non-Hispanic Black' then order=_N_-1;
    if by='Non-Hispanic White' then order=_N_-3;
    if by='Other/Unknown' then delete;
	if by='<=180 days' then order=_N_-1.5; 
	if by='6-17 yo' then do; 
		order=_N_-2.5; Exp_IDR='N/A'; Unexp_IDR='N/A'; 
		UnadjRR='N/A'; AdjRR='N/A'; UnadjVE='N/A'; AdjVE='N/A'; end;
    drop adjVE_adjCI separator4; *dropped adjVE_adjCI for this study as no multi-testing adjustment is planned;
    run;
    proc sort data=table3; by order; run;

	title 'Table 3';
    proc print data=table3; run;
	title;

***************************************************************;
** Table 4. rVE against medically attended SARS-CoV-2 infection and COVID-19 hospital death;

* Medically attended infection (Not adjusted for antiviral);
	%let varlist=age_grp sex race_eth index_month va_grp preventive_care ndose_c days_since_last_vac_c medical_center_and_area hist_covid_diag; 
    %cohort_analysis(data_analysis, covid_infec_all, &exp., &unexp., py_infec_all, lpy_infec_all, py_infec_all_mons, infec_all_rVE, survivalplot='Y');
    %cohort_analysis(data_analysis, covid_infec_edu, &exp., &unexp., py_infec_edu, lpy_infec_edu, py_infec_edu_mons, infec_edu_rVE, survivalplot='Y');
* Hosp death (Collapsed age_grp. Dropped med center);
	%let varlist=age_grp_comb3 sex race_eth index_month va_grp preventive_care ndose_c days_since_last_vac_c antiviral hist_covid_diag; 
	%cohort_analysis(data_analysis, covid_hosp_death, &exp., &unexp., py_dead, lpy_dead, py_dead_mons, dead_rVE, survivalplot='Y');

* Set results together;
    data blank; by=' '; run;
    data table4; 
    set Rslt_infec_all_rVE Rslt_infec_edu_rVE
		Rslt_dead_rVE;
    drop adjVE_adjCI separator4; *dropped adjVE_adjCI for this study as no multi-testing adjustment is planned;
    run;

	title 'Table 4';
    proc print data=table4; run;
	title;
