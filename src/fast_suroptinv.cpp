#include <RcppArmadillo.h>
// #include <iostream>
// using namespace std;
using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec phi(arma::vec myarg)
{
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
    double invr2 = 1/sqrt(2);
    double x;
    //myarg is a column vector

    long n = myarg.n_rows;
    arma::vec result(n);

    for(long ii=0;ii<n;ii++){
        // Save the sign of x
        x = myarg(ii);

        int sign = 1;
        if (x < 0)
            sign = -1;

        x = fabs(x)*invr2; ///sqrt(2.0);

        // A&S formula 7.1.26
        double t = 1.0/(1.0 + p*x);
        //double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
        double y = 1.0 - (a1 + a2*t + a3*t*t + a4*t*t*t + a5*t*t*t*t) *t*exp(-x*x);

        result(ii) = 0.5*(1.0 + sign*y);
    }

    return(result);
}


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::vec fast_suroptinv(
        const long N, const long p,const long nsimu,
        const arma::mat& allsimu,
        const arma::mat& allsimu_centered,
        const arma::mat& kn,
        const arma::mat& allKninv,
        const arma::mat& lambda,
        const double mn,const double sn,const double T,
        const arma::mat& randmatrix
) {

    // allsimu: p x nsimu*N          => N batches of nsimu simulations in p points
    // allsimu_centered: p x nsimu*N => N batches of nsimu simulations in p points
    // kn: p x N                     => a set of N vectors of size p : kriging covariances between xn+1 and the simu points
    // allKn.inv: p*N x p            => N covariances matrices p x p
    // lambda: p x N                 => a set of N vectors of kriging weights of xn+1 for prediction in p points
    // randmatrix: N x nsimu         => some randoms N(0,1)

    //declarations
    long index0, index, nvalid;
    double current_ui, current_vi, ui_sort_ind1, vi_sort_ind2, infinity, minus_infinity;
    infinity = std::numeric_limits<double>::max();
    minus_infinity = -1*infinity;

    arma::mat tkn_Kninv(N,p);
    arma::colvec fact1(p);
    arma::rowvec fact2(p);
    arma::colvec lambda_i(p);

    arma::mat m_nplusp(N,nsimu);
    arma::colvec s_nplusp(N);
    arma::rowvec z_xnplus1(nsimu);
    arma::rowvec v(nsimu);

    arma::vec ai(p);
    arma::vec ui(nsimu);
    arma::vec vi(nsimu);
    arma::vec validuivi(nsimu);
    arma::vec big_prob(N);

    // int start = clock();

    // first multiply the i-th column of kn by the matrix number i in allKninv
    for(long i=0; i<N; i++) {
        fact1 = kn.col(i);
        index0 = i*p;
        for(long j=0; j<p; j++) {
            index = index0 + j;
            fact2 = allKninv.row(index);
            tkn_Kninv(i,j) = arma::sum(fact2*fact1);
            // note * is the matrix product. Using * with two row (or column) vectors does not work
        }
    }

    // int end = clock();
    // cout << "100    ";
    // cout << end-start << endl;
    // start = clock();

    // then compute updated kriging means and variances to simulate conditionaly
    // on n+p points
    // the updated kriging mean m_nplusp of xnew depends on the simulation j
    // the updated kriging variance : no
    double sn2 = sn*sn;
    for(long i=0; i<N; i++) {
        fact2 = tkn_Kninv.row(i);
        for(long j=0;j<nsimu;j++){
            index=i*nsimu+j;
            m_nplusp(i,j) = mn+arma::sum(fact2*allsimu_centered.col(index));
        }
        s_nplusp(i) = sn2 - arma::sum(fact2*kn.col(i));
    }

    // end = clock();
    // cout << "117    ";
    // cout << (end-start) << endl;
    // start = clock();

    // then simulate and compute ai, ui,vi, and bigprob[i]
    // this loop goes until the end of the function
    for(long i=0; i<N; i++) {

        // deal with the simulations of integration point i
        v = randmatrix.row(i);
        z_xnplus1 = m_nplusp.row(i) + sqrt(s_nplusp(i)) * v;
        lambda_i = lambda.col(i);
        arma::colvec invlambda_i = 1.0 / lambda_i;

        ui = minus_infinity + arma::zeros(nsimu);
        vi = infinity + arma::zeros(nsimu);
        validuivi = arma::zeros(nsimu);

        index0=i*nsimu;
        for(long j=0; j<nsimu; j++) {
            index=index0+j;
            ai = (T - allsimu.col(index)) % invlambda_i + z_xnplus1(j);
            // This is the threshold called u_i or v_i in Chevalier's PhD manuscript.
            // which corresponds (in the manuscript) to the formula (T - a_i^j) / b^j
            // The name ai for this variable is somehow missleading.

            current_vi = infinity;
            current_ui = minus_infinity;

            for(long pp=0; pp<p;pp++){
                double lambdapp = lambda_i(pp);
                if(lambdapp>0){
                    current_vi = std::min(current_vi,ai(pp)); // upper bound vi
                }else{
                    current_ui = std::max(current_ui,ai(pp)); // lower bound vi
                }
            }

            ui(j) = current_ui;
            vi(j) = current_vi;
            validuivi(j) = (current_vi>current_ui); // valid if the inteval exists
        }

        // end = clock();
        // cout << "192        ";
        // cout << end-start << endl;
        // start = clock();

        nvalid = arma::sum(validuivi); // number of simu where we may not exceed T
        arma::vec validui(nvalid);
        arma::vec validvi(nvalid);
        long count=0;

        // copy the intervals [ui,vi] when they exist
        for(long pp=0;pp<nsimu;pp++){
            if( (vi(pp) > ui(pp)) ){
                validui(count) = ui(pp);
                validvi(count) = vi(pp);
                count++;
            }
        }

        // at this time we still did not normalize by substracting mn(xn+1)
        // and dividing by sn(xn+1)

        arma::vec ui_sort(nvalid);
        arma::vec vi_sort(nvalid);

        ui_sort = arma::sort(validui);
        vi_sort = arma::sort(validvi);

        // end = clock();
        // cout << "220        ";
        // cout << end-start << endl;
        // start = clock();

        arma::vec ui_and_vi(2*nvalid);
        arma::vec is_ui(2*nvalid);
        long ind1=0;
        long ind2=0;

        // merges ui and vi in a single sorted array
        // by recording if the number is an ui or not
        for(long uu=0;uu<(2*nvalid);uu++){
            if(ind1>=nvalid){
                ui_sort_ind1 = infinity;
            }else{
                ui_sort_ind1 = ui_sort(ind1);
            }
            if(ind2>=nvalid){
                vi_sort_ind2 = minus_infinity;
            }else{
                vi_sort_ind2 = vi_sort(ind2);
            }

            if( ui_sort_ind1<vi_sort_ind2 ){
                // record ui_sort_ind1 and increment ind1
                is_ui(uu) = 1;
                ui_and_vi(uu) = ui_sort_ind1;
                ind1++;
            }else{
                is_ui(uu) = -1;
                ui_and_vi(uu) = vi_sort_ind2;
                ind2++;
            }
        }

        // end = clock();
        // cout << "251        ";
        // cout << end-start << endl;
        // start = clock();

        arma::vec prop_simu(2*nvalid);
        arma::vec pn_oneminuspn(2*nvalid);
        arma::vec pnorms(2*nvalid);

        // This gives the proportion of simulation that don't exceed T, in function of the response
        prop_simu = arma::cumsum(is_ui)/nsimu;
        pn_oneminuspn = prop_simu % (1-prop_simu); // % does the same as * in R. Element per element product

        // end = clock();
        // cout << "270        ";
        // cout << end-start << endl;
        // start = clock();

        ui_and_vi = (ui_and_vi - mn)/sn;
        pnorms = phi(ui_and_vi);

        // end = clock();
        // cout << "279        ";
        // cout << end-start << endl;
        // start = clock();

        for(long uu=0;uu<(2*nvalid-1);uu++){
            // the prop_simu vector re-used here
            prop_simu(uu) = pn_oneminuspn(uu) * ( pnorms(uu+1) -  pnorms(uu)  );
        }

        big_prob(i) = arma::sum(prop_simu);
    }

    // end = clock();
    // cout << "304    ";
    // cout << end-start << endl;
    // start = clock();

    return( big_prob );
}

