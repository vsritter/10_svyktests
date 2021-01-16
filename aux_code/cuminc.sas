
* Calculation of cumulative incidence curves, associated standard error, pointwise
* confidence limits, and Gray test if comparing groups
*
* Usage:
*   %cuminc(
*        evtvar   = ,       Event variable - must be numeric, required
*                           Each value represents a distinct event, for example,
*                           0 = Censored (defined as such below), 1 = Event #1, 2 = Event #2,
*                           etc.  No limit on the number of competing events.
*                           Note that events occurring at time 0 are included in the calculation.
*                           If the user wishes to exclude these events, they should not appear
*                           in the dataset passed to the macro.
*        timevar  = ,       Survival time variable - required
*        grpvar   = ,       By-group variable - must be numeric if specified, 
*                           default: No group (All records)
*                           If specified, group-specific estimates are produced and the Gray
*                           test for group differences is performed.
*        stratvar = ,       Stratification variable, only applicable when a group variable has
*                           been specified.  Allows the user to perform a stratified group test.
*                           Must be numeric if specified, default: No strata 
*                           Note: No strata-specific estimates are produced
*        cencode  = ,       Code for censored observations - limited to a single value, default: 0
*        alpha    = ,       Alpha to be used for confidence limit calculation, default: 0.05
*        rho      = ,       Rho value for Gray test, specifies form of difference to be tested
*                           0 = Hazard ratio, 1 = Odds ratio, 2 = Cumulative risk ratio
*                           See Gray, 1988 for more detail.  Default: 0
*        inds     = ,       Input dataset - default: _LAST_
*        outds    = ,       Output dataset for estimates + variance and confidence limits 
*                           default: _OUTCI
*        statds   = ,       Output dataset for Gray test statistics and p-values - default: _OUTSTAT
*        timelist = ,       List of times for which to display CIs (as: t1 t2 t3 t4)
*        reduceout= ,       (Y|N) Limit &outds to only times in timelist - default: N
*        printit  = ,       (Y|N) Display estimates (can set to N if only output ds required) -
*                           default: Y
*   )
*
* Date: 25Mar2012
* Rev:  28Mar2012  Report SE instead of variance
*                  Fixed divide-by-0 error in CI calculation when CumInc = 0
*       08Nov2012  Fixed incorrect handling of censor codes when specified to not equal 0
*                  Added "requested time" column when timelist is specified
*       29Nov2012  Bypass CI construction if no events occurred
*       30May2013  Misc documentation changes
* 
* Author: 
*    Brad Hammill
*    brad.hammill@duke.edu
*    [adapted from R code for cmprsk.cuminc function]
* 
* References for cumulative incidence function and Gray test, if needed:
*    Kalbfleisch JD, Prentice RL. (1980). The Statistical Analysis of Failure Time Data. New York: 
*    John Wiley & Sons, Inc.
*
*    Gray RJ. (1988). A Class of K-Sample Tests for Comparing the Cumulative Incidence of
*    a Competing Risk.  Annals of Statistics, 16(3):1141-54.
;


%macro cuminc(
    evtvar   = ,
    timevar  = ,
    grpvar   = _allpts_,
    stratvar = _strata_,
    cencode  = 0,
    rho      = 0,
    alpha    = .05,
    inds     = _LAST_,
    outds    = _outci,
    statds   = _outstat,
    timelist = -9,
    reduceout = N,
    printit  = Y

);

    data _ci_;
        set &inds;

        _allpts_ = 1;
        _strata_ = 1;
    run;

    proc sort data=_ci_;
        by &timevar;
    run;

    proc iml;

        * Gray test within strata mechanics;

        start crst(y, m, ig, n, ng, rho, s, v, ng1);
            
            a   = J(ng, ng, 0);
            skm = J(ng, 1, 1);
            skmm= J(ng, 1, 1);
            f1m = J(ng, 1, 0);
            f1  = J(ng, 1, 0);
            rs  = J(ng, 1, 0);  
            v2  = J(ng1, ng, 0);
            v3  = J(ng, 1, 0);
            c   = J(ng, ng, 0);

            * Initialize rs w/# obs per group;
            do i = 1 to ng;
                rs[i] = ncol(loc(ig = i));
            end;

            fm = 0;
            f = 0;
            * begin looping over unique event times;
            doidx = loc(y[1:nrow(y) - 1] ^= y[2:nrow(y)]) || nrow(y);
            ll = 1;
            do lu = 1 to ncol(doidx);
                * d will contain the # in each group censored (=1), failed from
                * cause 1 (=2), and failing from other causes (=3);
                d = J(3, ng, 0);
                nd1 = 0;
                nd2 = 0;

                do i = ll to doidx[lu];
                    j = ig[i];
                    k = m[i] + 1;
                    d[k, j] = d[k, j] + 1;
                end;
                do i = 1 to ng;
                    nd1 = nd1 + d[2, i];
                    nd2 = nd2 + d[3, i];
                end;

                if sum(nd1, nd2) > 0 then do;
                    tr = 0;
                    tq = 0;
                    do i = 1 to ng;
                        if (rs[i] > 0) then do;
                            td = d[2, i] + d[3, i];
                            * skmm is left continuous,  and skm right continuous,  km est.;
                            skm[i] = skmm[i] * (rs[i] - td) / rs[i];
                            * f1m is left continuous,  and f1 right continuous,  cuminc est.;
                            f1[i] = f1m[i] + (skmm[i] * d[2, i]) / rs[i];
                            * in notation of the paper,  tr is \sum_r\hat{h}_r,  and tq is \sum_r R_r;
                            tr = tr + rs[i] / skmm[i];
                            tq = tq + rs[i] * (1 - f1m[i]) / skmm[i];
                        end;
                    end;
                    f = fm + nd1 / tr;
                    fb = (1 - fm)**rho;
                    do i = 1 to ng;
                        do j = i to ng;
                            a[i, j] = 0;
                        end;
                        if (rs[i] > 0) then do; 
                            t1 = rs[i] / skmm[i];
                            a[i, i] = fb * t1 * (1 - t1 / tr);
                            c[i, i] = c[i, i] + a[i, i] * nd1 / (tr * (1 - fm));
                            k = i + 1;
                            if (k <= ng) then do;
                                do j = k to ng;
                                    if (rs[j] > 0) then do;
                                        a[i, j] =  - fb * t1 * rs[j] / (skmm[j] * tr);
                                        c[i, j] = c[i, j] + a[i, j] * nd1 / (tr * (1 - fm));
                                    end;
                                end;
                            end;
                        end;
                    end;
                    do i = 2 to ng;
                        k = i - 1;
                        do j = 1 to k;
                            a[i, j] = a[j, i];
                            c[i, j] = c[j, i];
                        end;
                    end;
                    do i = 1 to ng1;
                        if (rs[i] > 0) then do;
                            s[i] = s[i] + fb * (d[2, i] - nd1 * rs[i] * (1 - f1m[i]) / (skmm[i] * tq));
                        end;
                    end;
                    if (nd1 > 0) then do;
                        do k = 1 to ng;
                            if (rs[k] > 0) then do;
                                t4 = 1;
                                if (skm[k] > 0) then t4 = 1 - (1 - f) / skm[k];
                                t5 = 1;
                                if (nd1 > 1) then t5 = 1 - (nd1 - 1) / (tr * skmm[k] - 1);
                                t3 = t5 * skmm[k] * nd1 / (tr * rs[k]);
                                v3[k] = v3[k] + t4 * t4 * t3;
                                do i = 1 to ng1;
                                    t1 = a[i, k] - t4 * c[i, k];
                                    v2[i, k] = v2[i, k] + t1 * t4 * t3;
                                    do j = 1 to i;
                                        l = i * (i - 1) / 2 + j;
                                        t2 = a[j, k] - t4 * c[j, k];
                                        v[l] = v[l] + t1 * t2 * t3;
                                    end;
                                end;
                            end;
                        end;
                    end;
                    if (nd2 > 0) then do;
                        do k = 1 to ng;
                            if (skm[k] > 0 & d[3, k] > 0) then do;
                                t4 = (1 - f) / skm[k];
                                t5 = 1;
                                if (d[3, k] > 1) then t5 = 1 - (d[3, k] - 1.0) / (rs[k] - 1.0);
                                t3 = t5 * ((skmm[k] ** 2) * d[3, k]) / (rs[k] ** 2);
                                v3[k] = v3[k] + t4 * t4 * t3;
                                do i = 1 to ng1;
                                    t1 = t4 * c[i, k];
                                    v2[i, k] = v2[i, k] - t1 * t4 * t3;
                                    do j = 1 to i;
                                        l = i * (i - 1) / 2 + j;
                                        t2 = t4 * c[j, k];
                                        v[l] = v[l] + t1 * t2 * t3;
                                    end;
                                end;
                            end;
                        end;
                    end;
                end;
                do i = ll to doidx[lu];
                    j = ig[i];
                    rs[j] = rs[j] - 1;
                end;
                fm = f;
                do i = 1 to ng;
                    f1m[i] = f1[i];
                    skmm[i] = skm[i];
                end;
                ll = doidx[lu] + 1;
            end;
            l = 0;
            do i = 1 to ng1;
                do j = 1 to i;
                    l = l + 1;
                    do k = 1 to ng;
                        v[l] = v[l] + c[i, k] * c[j, k] * v3[k];
                        v[l] = v[l] + c[i, k] * v2[j, k];
                        v[l] = v[l] + c[j, k] * v2[i, k];
                    end;
                end;
            end;
        
            return;

        finish; /* crst */


        **** START MAIN ****;
    
        use _ci_;

        * Read in all evt, group, timing, strata data;
        read all var {&evtvar} into evt;
        read all var {&grpvar} into grp;
        read all var {&timevar} into tm;
        read all var {&stratvar} into strat;


        * Counts of things;
        nobs = nrow(evt);

        * Current and any event vector -- No need to assign censor code to 0 since init = 0;
        cur = j(nobs, 1, 0);
        any = j(nobs, 1, 0);

        * Unique non-censor event value, group values, and times;
        evtlev = setdif(unique(evt), &cencode);
        grplev = unique(grp);
        tmlev  = unique(tm);
        strlev = unique(strat);

        ngrp = ncol(grplev);
        nevt = ncol(evtlev);
        ntm  = ncol(tmlev);
        nstr = ncol(strlev);

        * Recode group;
        _grp = J(nobs, 1, .);
        do i = 1 to ngrp;
            _grp[loc(grp = grplev[i])] = i;
        end;

        * Make timelist;
        tl = { &timelist };

        * Assign 1 to any event in any;
        evtidx = loc(evt ^= &cencode);
        any[evtidx] = 1;

        * Step through events;
        do e = 1 to nevt; 

            * Define current event vector;
            evtidx = loc(evt = evtlev[e]);
            cur[evtidx] = 1;

            if loc(evt ^= evtlev[e] & any = 1) then do;
                othidx = loc(evt ^= evtlev[e] & any = 1);
                cur[othidx] = 2;
            end;

            **** BEGIN ACTUAL GRAY TEST FOR EVENT E ****;

            if ngrp > 1 then do;

                * Calculate score and variance matrix for comparing CI curves;
                    
                ng1 = ngrp - 1;
                ng2 = ngrp * ng1 / 2;
                ng3 = 4 * ngrp + 1;
                ng4 = ngrp * ngrp;
                ng5 = 4 * ngrp;
                ng6 = ngrp * (4 + 3 * ngrp);
                rho = &rho ;

                v   = J(ng2, 1, 0);
                s   = J(ng1, 1, 0);
                vs  = J(ng1, ng1, 0);

                do ks = 1 to nstr;
                    stridx = loc(strat = strlev[ks]);
                    nks = ncol(stridx);

                    * Strata-specific tm, evt, grp vectors;
                    ys = tm[stridx];
                    ms = cur[stridx];
                    igs = _grp[stridx];

                    * Re-init vt and st;
                    vt  = J(ng2, 1, 0);
                    st  = J(ng1, 1, 0);

                    run crst(ys, ms, igs, nks, ngrp, rho, st, vt, ng1);

                    l = 0;
                    do i = 1 to ng1;
                        s[i] = s[i] + st[i];
                        do j = 1 to i;
                            l = l + 1;
                            v[l] = v[l] + vt[l];
                        end;
                    end;
                end;
                l = 0;
                do i = 1 to ng1;
                    do j = 1 to i;
                        l = l + 1;
                        vs[i, j] = v[l];
                        vs[j, i] = vs[i, j];
                    end;
                end;

                _q = s` * inv(vs) * s;
                _p = 1 - probchi(_q, ng1);

                _stat = evtlev[e] || _q || _p || ng1;

                * Append results;
                if e = 1 then
                    statout = _stat;
                else
                    statout = statout // _stat;
            end;

            **** END ACTUAL GRAY TEST FOR EVENT I ****;


            * Step through groups;
            do g = 1 to ngrp;

                * Get group-specific values;
                grpidx = loc(grp = grplev[g]);
                atrisk = ncol(grpidx);

                _cur = cur[grpidx];
                _tm = tm[grpidx];
                _tmlev = unique(_tm);
                _ntm = ncol(_tmlev); *******;

                * Create null output vectors;
                _out = j(_ntm, 9, .);
                _out[,1] = grplev[g]; * Group value;
                _out[,2] = evtlev[e]; * Event value;
                _out[,3] = _tmlev`;      * Event times;
                * 4 = CumInc, 5 = SE, 6 = LCL, 7 = UCL, 8 = ATRISK, 9 = EVENTS;

                * Figure rows to show and append to _out;
                if tl = -9 then do;
                    show = J(_ntm, 1, 1);
                    _out = _out || show;
                end;
                else do;
                    show = J(_ntm, 1, 0);
                    showtime = J(_ntm, 1, .);

                    do i = 1 to ncol(tl);
                        show[max(loc(_out[,3] <= tl[i]))] = 1;
                        showtime[max(loc(_out[,3] <= tl[i]))] = tl[i];
                    end;

                    _out = _out || show;
                    _out = _out || showtime;
                end;


                * For later pre-pending;
                if tl = -9 then 
                    _out0 = j(1, 10, 0); 
                else 
                    _out0 = j(1, 11, 0); * Requires an addl column for the requested time labels;

                _out0[1] = grplev[g];
                _out0[2] = evtlev[e];
                _out0[8] = atrisk;


                **** BEGIN GRAY MACRO ADAPTATION ****;

                * Specific edits/translations from original Gray macro:
                *   y   > _tm
                *   ic  > _any
                *   icc > _cur
                *   n   > atrisk
                ;

                * Identify last lines for unique times;
                doidx = loc(_tm[1:nrow(_tm) - 1] ^= _tm[2:nrow(_tm)]) || nrow(_tm);

                fk = 1;
                nf = 0;
                v1 = 0;
                v2 = 0;
                v3 = 0;
                f = 0;
                v = 0;
                rs = atrisk;
                ll = 1;

                do l = 1 to ncol(doidx);
                    nd1 = ncol(loc(_cur[ll:doidx[l]] = 1)); 
                    nd2 = ncol(loc(_cur[ll:doidx[l]] = 2)); 

                    nd = nd1 + nd2;
                  
                    if (nd ^= 0) then do;
                        fkn = fk * (rs - nd) / rs;

                        if (nd1 > 0) then do;
                            f = f + fk * nd1 / rs;
                        end;
                
                        if nd2 > 0 then if fkn > 0 then do;
                            t5 = 1;
                            if (nd2 > 1) then t5 = 1 - (nd2 - 1.) / (rs - 1.);
                            t6 = fk * fk * t5 * nd2 / (rs * rs);
                            t3 = 1. / fkn;
                            t4 = f / fkn;
                            v1 = v1 + t4 * t4 * t6;
                            v2 = v2 + t3 * t4 * t6;
                            v3 = v3 + t3 * t3 * t6;
                        end;

                        if (nd1 > 0) then do;
                            t5 = 1;
                            if (nd1 > 1) then t5 = 1 - (nd1 - 1.) / (rs - 1.);
                            t6 = fk * fk * t5 * nd1 / (rs * rs);
                            t3 = 0;
                            if (fkn > 0) then t3 = 1 / fkn;
                            t4 = 1 + t3 * f;
                            v1 = v1 + t4 * t4 * t6;
                            v2 = v2 + t3 * t4 * t6;
                            v3 = v3 + t3 * t3 * t6;
                            t2 = f;
                            v = v1 + t2 * t2 * v3 - 2 * t2 * v2;
                        end; 
                        
                        fk = fkn;
                        nf = nf + nd1;
                    end; 

                    _out[l, 4] = f;
                    _out[l, 5] = sqrt(v); 
                    
                    _out[l, 8] = rs;
                    _out[l, 9] = nf;

                    rs = atrisk - doidx[l];
                    ll = doidx[l] + 1;
                end;

                * Calculate UCL/LCL;
                vidx = loc(_out[, 4] > 0);
                if nrow(vidx) > 0 then do;
                    za = probit(1 - (&alpha / 2));
                    _out[vidx, 6] = _out[vidx, 4] # exp(-1 # za # _out[vidx, 5] / _out[vidx, 4]);
                    _out[vidx, 7] = _out[vidx, 4] # exp(+1 # za # _out[vidx, 5] / _out[vidx, 4]);
                end;

                **** END GRAY MACRO ADAPTATION ****;

                * Prepend t = 0 info results;
                _out = _out0 // _out;

                * Append results;
                if g = 1 & i = 1 then
                    ciout = _out;
                else
                    ciout = ciout // _out;
            end;

        end;

        if ngrp > 1 then do;
            * Output statistics;
            statlbl = {"EVENT" "STAT" "P" "DF"};
            create &statds from statout [ colname = statlbl ];
            append from statout;
        end;

        * Output estimates;
        if tl = -9 then 
            estlbl = {"GROUP", "EVENT", "TIME", "CUMINC", "STDERR", "LCL", "UCL", "ATRISK", "EVENTS", "SHOW"};
        else 
            estlbl = {"GROUP", "EVENT", "TIME", "CUMINC", "STDERR", "LCL", "UCL", "ATRISK", "EVENTS", "SHOW", "REQTIME"};

        create &outds from ciout [ colname = estlbl ];

        * If REDUCEOUT, then only include where show = 1;
        %if %upcase(&reduceout) = Y %then %do;
            cilimit = ciout[loc(ciout[,10] = 1),];
            append from cilimit;
        %end;
        %else %do;
            append from ciout;
        %end;
    quit;

    * Print stuff;
    %if &printit = Y %then %do;
        %if %sysfunc(exist(&statds)) = 1 %then %do;
            proc print noobs data=&statds;
                format p pvalue8.4;
            run;
        %end;

        proc print noobs data=&outds;
            by event group;
            id event group;
            where show;
            var TIME CUMINC STDERR LCL UCL ATRISK EVENTS;
            format 
                cuminc lcl ucl 8.6 stderr 9.7
            ;
        run;
    %end;

    * Delete the show indicator (?);
    data &outds;
        set &outds;
    run;

    * Delete the single non-returned dataset;
    proc delete data=_ci_;
    run;

%mend;



