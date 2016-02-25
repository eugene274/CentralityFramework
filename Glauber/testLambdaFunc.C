void testLambdaFunc() { 


   TH1D * h1 = new TH1D("h1","h1",100,-5,5); 
   h1->FillRandom("gaus"); 

   TGraph * g = new TGraph(h1); 

   TF1 * f1 = new TF1("f",[&](double *x, double *p) { return p[0]*g->Eval(x[0] ); },-5,5, 1); 

   h1->Fit(f1); 

   // this fixes the problem
   h1->Fit(f1,"N"); 

   TF1 * fittedFunc = (TF1*) f1->Clone(); 
   h1->GetListOfFunctions()->Add(fittedFunc); 

   h1->Draw(); 
}
