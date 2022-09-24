void WCSimVertexFinder::ConePropertiesLnL(Double_t coneParam0, Double_t coneParam1, Double_t coneParam2, Double_t& coneAngle, Double_t& coneFOM)
{  
  // nuisance parameters
  // ===================
  Double_t alpha  = coneParam0;
  Double_t alpha0 = coneParam1;
  Double_t beta   = coneParam2;

  // internal variables
  // ==================
  Double_t deltaAngle = 0.0;
  Double_t sigmaAngle = 7.0;
  Double_t deltaAngle0 = 42.0*alpha0;
  
  Double_t digitQ = 0.0;
  Double_t sigmaQmin = 1.0;
  Double_t sigmaQmax = 10.0;
  Double_t sigmaQ = 0.0;

  Double_t A = 0.0;
  
  Double_t PconeA = 0.0;
  Double_t PconeB = 0.0;
  Double_t Pmu = 0.0;
  Double_t Pel = 0.0;

  Double_t Pcharge = 0.0;
  Double_t Pangle = 0.0;
  Double_t P = 0.0;

  Double_t chi2 = 0.0;
  Double_t ndof = 0.0;

  Double_t angle = 46.0;
  Double_t fom = 0.0;

  // hard-coded parameters: 200 kton (100 kton)
  // ==========================================
  Double_t lambdaMuShort = 0.5; //  0.5;
  Double_t lambdaMuLong  = 5.0; // 15.0;
  Double_t alphaMu =       1.0; //  4.5;

  Double_t lambdaElShort = 1.0; //  2.5;
  Double_t lambdaElLong =  7.5; // 15.0;
  Double_t alphaEl =       6.0; //  3.5;

  // numerical integrals
  // ===================
  fSconeA = 21.9938;  
  fSconeB =  0.0000;

  // inside cone
  Int_t nbinsInside = 420;
  for( Int_t n=0; n<nbinsInside; n++ ){
    deltaAngle = -42.0 + (n+0.5)*(42.0/(double)nbinsInside);
    fSconeB += 1.4944765*sin( (42.0+deltaAngle)*(TMath::Pi()/180.0) )
                           *( 1.0/(1.0+(deltaAngle*deltaAngle)/(deltaAngle0*deltaAngle0)) )
                           *( 42.0/(double)nbinsInside );
  }

  // outside cone
  if( fIntegralsDone == 0 ){
    fSmu = 0.0;
    fSel = 0.0;

    Int_t nbinsOutside = 1380;
    for( Int_t n=0; n<nbinsOutside; n++ ){
      deltaAngle = 0.0 + (n+0.5)*(138.0/(double)nbinsOutside);

      fSmu += 1.4944765*sin( (42.0+deltaAngle)*(TMath::Pi()/180.0) )
                          *( 1.0/(1.0+alphaMu*(lambdaMuShort/lambdaMuLong)) )*( 1.0/(1.0+(deltaAngle*deltaAngle)/(lambdaMuShort*lambdaMuShort)) 
                                            + alphaMu*(lambdaMuShort/lambdaMuLong)/(1.0+(deltaAngle*deltaAngle)/(lambdaMuLong*lambdaMuLong)) )
                          *( 138.0/(double)nbinsOutside );

      fSel += 1.4944765*sin( (42.0+deltaAngle)*(TMath::Pi()/180.0) )
                          *( 1.0/(1.0+alphaEl*(lambdaElShort/lambdaElLong)) )*( 1.0/(1.0+(deltaAngle*deltaAngle)/(lambdaElShort*lambdaElShort)) 
                                          + alphaEl*(lambdaElShort/lambdaElLong)/(1.0+(deltaAngle*deltaAngle)/(lambdaElLong*lambdaElLong)) )
                          *( 138.0/(double)nbinsOutside );
    }

    std::cout << " --- calculating integrals: Smu=" << fSmu << " Sel=" << fSel << std::endl;

    fIntegralsDone = 1;
  }


  // loop over digits
  // ================
  for( Int_t idigit=0; idigit<WCSimVertexGeometry::Instance()->GetNDigits(); idigit++ ){

    if( WCSimVertexGeometry::Instance()->IsFiltered(idigit) ){
      digitQ = WCSimVertexGeometry::Instance()->GetDigitQ(idigit);
      deltaAngle = WCSimVertexGeometry::Instance()->GetAngle(idigit) - 42.0;

      // pulse height distribution
      // =========================
      if( deltaAngle<=0 ){
        sigmaQ = sigmaQmax;
      }
      else{
        sigmaQ = sigmaQmin + (sigmaQmax-sigmaQmin)/(1.0+(deltaAngle*deltaAngle)/(sigmaAngle*sigmaAngle));
      }

      A = 1.0/(log(2.0)+0.5*TMath::Pi()*sigmaQ);

      if( digitQ<1.0 ){
        Pcharge = 2.0*A*digitQ/(1.0+digitQ*digitQ);
      }
      else{
        Pcharge = A/(1.0+(digitQ-1.0)*(digitQ-1.0)/(sigmaQ*sigmaQ));
      }

      // angular distribution
      // ====================
      A = 1.0/( alpha*fSconeA + (1.0-alpha)*fSconeB 
               + beta*fSmu + (1.0-beta)*fSel ); // numerical integrals 

      if( deltaAngle<=0 ){

        // pdfs inside cone:
        PconeA = 1.4944765*sin( (42.0+deltaAngle)*(TMath::Pi()/180.0) );
        PconeB = 1.4944765*sin( (42.0+deltaAngle)*(TMath::Pi()/180.0) )
                          *( 1.0/(1.0+(deltaAngle*deltaAngle)/(deltaAngle0*deltaAngle0)) );

        Pangle = A*( alpha*PconeA+(1.0-alpha)*PconeB );
      }         
      else{

        // pdfs outside cone
        Pmu = 1.4944765*sin( (42.0+deltaAngle)*(TMath::Pi()/180.0) )
                       *( 1.0/(1.0+alphaMu*(lambdaMuShort/lambdaMuLong)) )*( 1.0/(1.0+(deltaAngle*deltaAngle)/(lambdaMuShort*lambdaMuShort)) 
                                          + alphaMu*(lambdaMuShort/lambdaMuLong)/(1.0+(deltaAngle*deltaAngle)/(lambdaMuLong*lambdaMuLong)) );

        Pel = 1.4944765*sin( (42.0+deltaAngle)*(TMath::Pi()/180.0) )
                       *( 1.0/(1.0+alphaEl*(lambdaElShort/lambdaElLong)) )*( 1.0/(1.0+(deltaAngle*deltaAngle)/(lambdaElShort*lambdaElShort)) 
                                          + alphaEl*(lambdaElShort/lambdaElLong)/(1.0+(deltaAngle*deltaAngle)/(lambdaElLong*lambdaElLong)) );

        Pangle = A*( beta*Pmu+(1.0-beta)*Pel );
      }

      // overall probability
      // ===================
      P = Pcharge*Pangle;
      
      chi2 += -2.0*log(P);
      ndof += 1.0;
    }
  }

  // calculate figure of merit
  // =========================   
  if( ndof>0.0 ){
    fom = fBaseFOM - 5.0*chi2/ndof;
    angle = beta*43.0 + (1.0-beta)*49.0;
  }

  // return figure of merit
  // ======================
  coneAngle = angle;
  coneFOM = fom;

  return;
}