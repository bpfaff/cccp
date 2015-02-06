#include "cccp.h"
/*
 *
 * Methods for CONEC
 *
*/
using namespace arma;
/*
 * Inner product of two vectors in S.
*/
vec CONEC::sdot(mat s, mat z){
  vec ans = zeros<vec>(K);

  for(int i = 0; i < K; i++){
    if(cone[i] != "PSDC") {
      ans.at(i) = sdot_nlp(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
			   z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all));
    } else {
      ans.at(i) = sdot_s(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
			 z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
			 dims[i]);
    }
  }

  return ans;
}
/*
 * Norm of a vectors in S.
*/
double CONEC::snrm2(mat s){
  double ans = 0.0;

  for(int i = 0; i < K; i++){
    if(cone[i] != "PSDC"){
      ans += snrm2_nlp(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all));
    } else {
      ans += snrm2_s(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), dims[i]);
    } 
  }

  return ans;
}
/*
 * Product between two vectors in S.
*/
mat CONEC::sprd(mat s, mat z){
  mat ans(G.n_rows, 1);

  for(int i = 0; i < K; i++){
    if((cone[i] == "NLFC") || (cone[i] == "NNOC")){
      ans(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = 
	sprd_nl(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
		z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all));
    } else if(cone[i] == "SOCC"){
      ans(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = 
	sprd_p(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
	       z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all));
    } else if(cone[i] == "PSDC"){
      ans(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = 
	sprd_s(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
	       z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
	       dims[i]);
    }
  }

  return ans;
}
/*
 * One-element (neutral) with respect to a vector in S.
*/
mat CONEC::sone(){
  mat ans(G.n_rows, 1);

  for(int i = 0; i < K; i++){
    if((cone[i] == "NLFC") || (cone[i] == "NNOC")){
      ans(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = sone_nl(dims[i]);
    } else if(cone[i] == "SOCC"){
      ans(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = sone_p(dims[i]);
    } else if(cone[i] == "PSDC"){
      ans(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = sone_s(dims[i]);
    }
  }

  return ans;
}
/*
 * Inverse of product between two vectors in S.
*/
mat CONEC::sinv(mat s, mat z){
  mat ans(G.n_rows, 1);

  for(int i = 0; i < K; i++){
    if((cone[i] == "NLFC") || (cone[i] == "NNOC")){
      ans(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = 
	sinv_nl(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
		z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all));
    } else if(cone[i] == "SOCC"){
      ans(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = 
	sinv_p(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
	       z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all));
    } else if(cone[i] == "PSDC"){
      ans(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = 
	sinv_s(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
	       z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
	       dims[i]);
    }
  }

  return ans;
}
/*
 * Determining maximum step-size of a vector in S.
*/
vec CONEC::smss(mat u){
  vec ans = zeros<vec>(K);

  for(int i = 0; i < K; i++){
    if((cone[i] == "NLFC") || (cone[i] == "NNOC")){
      ans.at(i) = smss_nl(u(span(sidx.at(i, 0), sidx.at(i, 1)), span::all));
    } else if(cone[i] == "SOCC"){
      ans.at(i) = smss_p(u(span(sidx.at(i, 0), sidx.at(i, 1)), span::all));
    } else if(cone[i] == "PSDC"){
      ans.at(i) = smss_s(u(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), dims[i]);
    }
  }

  return ans;
}
/*
 * Applying maximum step-size to a vector in S (initial).
*/
mat CONEC::sams1(mat u, double alpha){
  mat temp;

  for(int i = 0; i < K; i++){
    temp = u(span(sidx.at(i, 0), sidx.at(i, 1)), span::all);
    if((cone[i] == "NLFC") || (cone[i] == "NNOC")){
      temp = sams1_nl(temp, alpha);
    } else if(cone[i] == "SOCC"){
      temp = sams1_p(temp, alpha);
    } else if(cone[i] == "PSDC"){
      temp = sams1_s(temp, alpha, dims[i]);
    }
    u(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = temp;
  }

  return u;
}
/*
Computation of Nesterov-Todd Scaling
*/
std::vector<std::map<std::string,mat> > CONEC::ntsc(mat s, mat z){
  std::vector<std::map<std::string,mat> > WList;
  std::map<std::string,mat> W;

  for(int i = 0; i < K; i++){
    if(cone[i] == "NLFC"){
      W = ntsc_n(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
		 z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all));
    } else if(cone[i] == "NNOC"){
      W = ntsc_l(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
		 z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all));
    } else if(cone[i] == "SOCC"){
      W = ntsc_p(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
		 z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all));
    } else if(cone[i] == "PSDC"){
      W = ntsc_s(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
		 z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), dims[i]);
    }
    WList.push_back(W);
    W.erase(W.begin(), W.end());
  }

  return WList;
}
/*
Update of Nesterov-Todd Scaling
*/
std::vector<std::map<std::string,mat> > CONEC::ntsu(mat s, mat z, std::vector<std::map<std::string,mat> > WList){ 
  std::map<std::string,mat> W;

  for(int i = 0; i < K; i++){
    if(cone[i] == "NLFC"){
      W = ntsu_n(WList[i], s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
		 z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all));
    } else if(cone[i] == "NNOC"){
      W = ntsu_l(WList[i], s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
		 z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all));
    } else if(cone[i] == "SOCC"){
      W = ntsu_p(WList[i], s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
		 z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all));
    } else if(cone[i] == "PSDC"){
      W = ntsu_s(WList[i], s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
		 z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), dims[i]);
    }
    WList[i] = W;
  }

  return WList;
}
/*
 * Scaling of vector in S by log-barrier function.
*/
mat CONEC::sslb(mat s, mat lambda, bool invers){
  mat ans(G.n_rows, 1);

  for(int i = 0; i < K; i++){
    if((cone[i] == "NLFC") || (cone[i] == "NNOC")){
      ans(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = 
	sslb_nl(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
		lambda(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), invers);
    } else if(cone[i] == "SOCC"){
      ans(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = 
	sslb_p(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
	       lambda(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), invers);
    } else if(cone[i] == "PSDC"){
      ans(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = 
	sslb_s(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
	       lambda(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), 
	       invers, dims[i]);
    }
  }

  return ans;
}
/*
 * Scaling of vector in S by Nesterov-Todd function.
*/
mat CONEC::ssnt(mat s, std::vector<std::map<std::string,mat> > WList, 
		bool invers, bool transp){

  for(int i = 0; i < K; i++){
    if(cone[i] == "NLFC"){
      s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = 
	ssnt_n(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), WList[i], invers);
    } else if(cone[i] == "NNOC"){
      s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = 
	ssnt_l(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), WList[i], invers);
    } else if(cone[i] == "SOCC"){
      s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = 
	ssnt_p(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), WList[i], invers);
    } else if(cone[i] == "PSDC"){
      s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = 
	ssnt_s(s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), WList[i], invers, 
	       transp);
    }
  }

  return s;
}
/*
Initial Nesterov-Todd scalings
*/
std::vector<std::map<std::string,mat> > CONEC::initnts(){
  std::vector<std::map<std::string,mat> > WList;
  std::map<std::string,mat> W;
  mat ans;

  for(int i = 0; i < K; i++){
    if(cone[i] == "NLFC"){
      ans = ones(dims[i],1);
      W.insert(std::pair<std::string,mat>("dnl", ans));
      W.insert(std::pair<std::string,mat>("dnli", ans));
      ans = zeros(dims[i],1);
      W.insert(std::pair<std::string,mat>("lambda", ans));
    } else if(cone[i] == "NNOC"){
      ans = ones(dims[i],1);
      W.insert(std::pair<std::string,mat>("d", ans));
      W.insert(std::pair<std::string,mat>("di", ans));
      ans = zeros(dims[i],1);
      W.insert(std::pair<std::string,mat>("lambda", ans));
    } else if(cone[i] == "SOCC"){
      ans = ones(1,1);
      W.insert(std::pair<std::string,mat>("beta", ans));
      ans = zeros(dims[i],1);
      ans.at(0,0) = 1.0;
      W.insert(std::pair<std::string,mat>("v", ans));
      ans = zeros(dims[i],1);
      W.insert(std::pair<std::string,mat>("lambda", ans));
    } else if(cone[i] == "PSDC"){
      ans = eye(dims[i],dims[i]);
      W.insert(std::pair<std::string,mat>("r", ans));
      W.insert(std::pair<std::string,mat>("rti", ans));
      ans = zeros(dims[i] * dims[i], 1);
      W.insert(std::pair<std::string,mat>("lambda", ans));
    }
    WList.push_back(W);
    W.erase(W.begin(), W.end());
  }

  return WList;
}
/*
Extracting Lagrange-Multipliers as matrix
*/
mat CONEC::getLambda(std::vector<std::map<std::string,mat> > WList){
  mat ans(G.n_rows, 1);

  for(int i = 0; i < K; i++){
    ans(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = WList[i]["lambda"]; 
  }

  return ans;
}
/*
Computation of: Sum of G_i'W_i^-1W_i^-1'G_i for i = 1, ..., K
*/
mat CONEC::gwwg(std::vector<std::map<std::string,mat> > WList){
  int n = G.n_cols;
  mat gwwg(n,n), temp(n,n), witg, wiwitg;
  gwwg.zeros();
  temp.zeros();

  for(int i = 0; i < K; i++){
    if(cone[i] == "NLFC"){
      witg = ssnt_n(G(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), WList[i], true);
      wiwitg = ssnt_n(witg, WList[i], true);
      temp = G(span(sidx.at(i, 0), sidx.at(i, 1)), span::all).t() * wiwitg;
    } else if(cone[i] == "NNOC"){
      witg = ssnt_l(G(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), WList[i], true);
      wiwitg = ssnt_l(witg, WList[i], true);
      temp = G(span(sidx.at(i, 0), sidx.at(i, 1)), span::all).t() * wiwitg;
    } else if(cone[i] == "SOCC"){
      witg = ssnt_p(G(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), WList[i], true);
      wiwitg = ssnt_p(witg, WList[i], true);
      temp = G(span(sidx.at(i, 0), sidx.at(i, 1)), span::all).t() * wiwitg;
    } else if(cone[i] == "PSDC"){
      witg = ssnt_s(G(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), WList[i], true, true);
      wiwitg = ssnt_s(witg, WList[i], true, false);
      temp = G(span(sidx.at(i, 0), sidx.at(i, 1)), span::all).t() * wiwitg;
    }
    gwwg = gwwg + temp;
  }

  return gwwg;
}
/*
Computation of: Sum of G_i'W_i^-1W_i^-1'G_i for i = 1, ..., K
*/
mat CONEC::gwwz(std::vector<std::map<std::string,mat> > WList, mat z){
  int n = G.n_cols;
  mat gwwz(n,1), temp(n,1), witz, wiwitz;
  gwwz.zeros();
  temp.zeros();

  for(int i = 0; i < K; i++){
    if(cone[i] == "NLFC"){
      witz = ssnt_n(z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), WList[i], true);
      wiwitz = ssnt_n(witz, WList[i], true);
      temp = G(span(sidx.at(i, 0), sidx.at(i, 1)), span::all).t() * wiwitz;
    } else if(cone[i] == "NNOC"){
      witz = ssnt_l(z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), WList[i], true);
      wiwitz = ssnt_l(witz, WList[i], true);
      temp = G(span(sidx.at(i, 0), sidx.at(i, 1)), span::all).t() * wiwitz;
    } else if(cone[i] == "SOCC"){
      witz = ssnt_p(z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), WList[i], true);
      wiwitz = ssnt_p(witz, WList[i], true);
      temp = G(span(sidx.at(i, 0), sidx.at(i, 1)), span::all).t() * wiwitz;
    } else if(cone[i] == "PSDC"){
      witz = ssnt_s(z(span(sidx.at(i, 0), sidx.at(i, 1)), span::all), WList[i], true, true);
      wiwitz = ssnt_s(witz, WList[i], true, false);
      temp = G(span(sidx.at(i, 0), sidx.at(i, 1)), span::all).t() * wiwitz;
    }
    gwwz = gwwz + temp;
  }

  return gwwz;
}
/*
Initializing PDV
*/
PDV* CONEC::initpdv(int p){
  PDV* pdv = new PDV();
  mat s(G.n_rows, 1);
  mat ans;

  pdv->x = zeros(n,1);
  pdv->y = zeros(p,1);
  for(int i = 0; i < K; i++){
    if((cone[i] == "NLFC") || (cone[i] == "NNOC")){
      ans = ones(dims[i], 1);
    } else if(cone[i] == "SOCC"){
      ans = zeros(dims[i], 1);
      ans.at(0,0) = 1.0;
    } else if(cone[i] == "PSDC") {
      ans = eye(dims[i],dims[i]);
      ans.reshape(dims[i] * dims[i], 1);
    } else {
      ans = zeros(dims[i], 1);
    }
    s(span(sidx.at(i, 0), sidx.at(i, 1)), span::all) = ans;
  }
  pdv->s = s;
  pdv->z = s;
  pdv->tau = 1.0;
  pdv->kappa = 1.0;

  return pdv;
}
/*
Updating Slack-variables
*/
mat CONEC::SorZupdate(mat SorZ, mat Lambda, double step){
  vec eval;
  mat tmpmat, evec;

  for(int j = 0; j < K; j++){
    if((cone[j] == "NLFC") || (cone[j] == "NNOC")){
      SorZ(span(sidx.at(j, 0), sidx.at(j, 1)), span::all) = \
	sams2_nl(SorZ(span(sidx.at(j, 0), sidx.at(j, 1)), span::all), step);
      SorZ(span(sidx.at(j, 0), sidx.at(j, 1)), span::all) = \
	sslb_nl(SorZ(span(sidx.at(j, 0), sidx.at(j, 1)), span::all), \
		Lambda(span(sidx.at(j, 0), sidx.at(j, 1)), span::all), true);
    } else if(cone[j] == "SOCC"){
      SorZ(span(sidx.at(j, 0), sidx.at(j, 1)), span::all) = \
	sams2_p(SorZ(span(sidx.at(j, 0), sidx.at(j, 1)), span::all), step);
      SorZ(span(sidx.at(j, 0), sidx.at(j, 1)), span::all) = \
	sslb_p(SorZ(span(sidx.at(j, 0), sidx.at(j, 1)), span::all), \
	       Lambda(span(sidx.at(j, 0), sidx.at(j, 1)), span::all), true);
    } else if(cone[j] == "PSDC"){
      tmpmat = SorZ(span(sidx.at(j, 0), sidx.at(j, 1)), span::all);
      tmpmat.reshape(dims[j], dims[j]);
      eig_sym(eval, evec, tmpmat);
      evec.reshape(dims[j] * dims[j], 1);
      SorZ(span(sidx.at(j, 0), sidx.at(j, 1)), span::all) = \
	sslb_s(evec, Lambda(span(sidx.at(j, 0), sidx.at(j, 1)), span::all), \
	       true, dims[j]);
      SorZ(span(sidx.at(j, 0), sidx.at(j, 1)), span::all) = \
	sams2_s(SorZ(span(sidx.at(j, 0), sidx.at(j, 1)), span::all), \
		step, Lambda(span(sidx.at(j, 0), sidx.at(j, 1)), span::all), \
		eval, dims[j]);
    }
  }

  return SorZ;
} 
