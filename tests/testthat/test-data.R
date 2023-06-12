# function to generate random testing datasets
gen_sim_data = function(m,p.x,p.t=0,n=0)
{
  XT.sim = lhs::create_oalhs(m,p.x+p.t,T,F)
  X.sim = XT.sim[,1:p.x]
  if(p.t>0){
    T.sim = XT.sim[,(p.x+1):(p.x+p.t)]
  } else{
    T.sim = NULL
  }
  # d_y = 10
  Y.sim = matrix(runif(m*10),nrow=10,ncol=m)

  if(n>0){
    X.obs = lhs::create_oalhs(n,p.x,T,F)
    Y.obs = matrix(runif(n*10),ncol=n)
  } else{
    X.obs = NULL
    Y.obs = NULL
  }
  return(list(X.sim=X.sim,T.sim=T.sim,Y.sim=Y.sim,X.obs=X.obs,Y.obs=Y.obs))
}

test_that("transform_xt works", {
  # test with only 1 x and 1 t as well as multiple - in the past the matrix operations have had issues with only a single column

  # no t - sim only
  test_dat = gen_sim_data(100,1)
  expect_error(FlaGP:::transform_xt(test_dat$X.sim), NA)
  test_dat = gen_sim_data(100,2)
  expect_error(FlaGP:::transform_xt(test_dat$X.sim), NA)  # t - sim + obs
  test_dat = gen_sim_data(100,1,1,5)
  expect_error(FlaGP:::transform_xt(test_dat$X.sim,test_dat$T.sim,test_dat$X.obs), NA)
  test_dat = gen_sim_data(100,2,2,5)
  expect_error(FlaGP:::transform_xt(test_dat$X.sim,test_dat$T.sim,test_dat$X.obs), NA)
})

test_that("transform_y works", {
  test_dat = gen_sim_data(100,1,1,5)
  # function already checks that mean = 0 and sd = 1 so we won't check that, just checking that it doesn't break
  expect_error(FlaGP:::transform_y(test_dat$Y.sim,y.ind.sim=NULL,Y.obs=NULL,y.ind.obs=NULL,center=T,scale=F,scaletype='scalar'), NA)
  expect_error(FlaGP:::transform_y(test_dat$Y.sim,y.ind.sim=matrix(seq_len(10)),Y.obs=NULL,y.ind.obs=NULL,center=T,scale=T,scaletype='rowwise'), NA)
  expect_error(FlaGP:::transform_y(test_dat$Y.sim,y.ind.sim=NULL,Y.obs=test_dat$Y.obs,y.ind.obs=NULL,center=F,scale=T,scaletype='scalar'), NA)
  expect_error(FlaGP:::transform_y(test_dat$Y.sim,y.ind.sim=matrix(seq_len(10)),Y.obs=NULL,y.ind.obs=NULL,center=F,scale=F,scaletype='rowwise'), NA)
})

test_that("basis works", {
  test_dat = gen_sim_data(100,2)
  # test with n.pc
  expect_error(FlaGP:::get_basis(test_dat$Y.sim,n.pc=3), NA)
  expect_error(FlaGP:::get_obs_basis(FlaGP:::get_basis(test_dat$Y.sim,n.pc=3),test_dat$Y.obs))
  # test with pct.var
  expect_error(FlaGP:::get_basis(test_dat$Y.sim,pct.var=.95), NA)
  expect_error(FlaGP:::get_obs_basis(FlaGP:::get_basis(test_dat$Y.sim,pct.var=.95),test_dat$Y.obs))
})

test_that("mv_lengthscales works", {

})
test_that("get_SC_inputs works", {

})
