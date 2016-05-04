
	// stop if fmax<0 (skipping phase re-centralization and the computation of coh)
	if( fmax < 0 ) return NaN;
	// when zeropha==true, shift phase range to [-pi/2, 3pi/2)
	if(zeropha) Pha.Transform( [&](float& val){if(val<-pio2)val+=twopi;} );
	// find phase median
	auto sigpha = Pha.sig.get(), sigcoh = Coh.sig.get();
	int ifmax = Pha.Index(fmax, true);
	std::vector<float> phaV; phaV.reserve(ifmax/2);
	for(int i=0; i<ifmax; i++) {
		if( sigcoh[i]<cohmin ) continue;	// exclude low coh points
		if( zeropha ) {	// and exclude points far away from 0/pi when zeropha==true
			float pha = sigpha[i]>=pio2 ? sigpha[i]-M_PI : sigpha[i];	// [-pi/2,pi/2)
			if( fabs(pha)>pio4 ) continue;
		}
		phaV.push_back(sigpha[i]);
	}
	float median = 0;
	if( phaV.size() > 0 ) {	// find median
		int imed = phaV.size()/2;
		std::nth_element(phaV.begin(), phaV.begin()+imed, phaV.end());
		median = phaV[imed];
	}
	// center-phase = median	if zeropha == false
	//					 = 0 or pi	if zeropha == true
	float phac = zeropha ? (fabs(median)<pio2?0:M_PI) : (median);
	// correct 2pi => [phac-pi, phac+pi)
	Pha.Transform( [&](float& val){
			if(val<phac-M_PI) val+=twopi;
			else if(val>=phac+M_PI) val-=twopi;
			} );
	/*
	// count for points near 0 or pi (with coh>cohmin)
	int count0 = 0, countpi = 0;
	for(int i=0; i<Pha.Index(fmax, true); i++) {
	if( sigcoh[i]<cohmin ) continue;
	float phabs = fabs(sigpha[i]);
	if( phabs<pio4 ) count0++;
	else if(phabs>pi3o4) countpi++;
	}
	int nvalid = std::max(count0, countpi);
	// theoretically (for tilt and compliance on Z), pha should be either 0 or pi:
	float phac = count0>countpi ? 0. : M_PI;
	// correct 2pi
	if( phac > M_PI*0.5 ) Pha.Transform( [&](float& val){if(val<0.)val+=twopi;} );
	*/
	// compute average coh between fb_avg and fe_avg Hz
	double cohavg = 0.;
	int i, ib = Pha.Index(fb_avg);
	for(i=ib; i<Pha.Index(fe_avg<fmax?fe_avg:fmax); i++) {
		float phadiff = fabs(sigpha[i]-phac);
		if(phadiff > pio4) continue;
		cohavg += sigcoh[i];
	}
	cohavg /= (i-ib);
	// save phac in Pha header
	Pha.shd.user1 = phac;
	//Coh.shd.user2 = tbeff; Coh.shd.user3 = teeff;
	return cohavg;

/*
std::vector<int> StaSacs::RemoveTilts( float Eperl, float Eperu, float tseg ) {
	_tseg = tseg;
	DetectNoiseWindows(Eperl, Eperu);
	// a vector to store (a total of 4) nfcorrs (which indicates the level of coherence between each channel)
	// to be returned
	std::vector<int> nfcorrV(4);
	// compute transfer function on both horizontals
	SacRec Coh1, Adm1, Pha1;
	int ncorr1 = CalcTransferF(sacZ, sacZamV, sacZphV, sacH1, sacH1amV, sacH1phV, Coh1, Adm1, Pha1);
	nfcorrV[0] = ncorr1;
	SacRec Coh2, Adm2, Pha2;
	int ncorr2 = CalcTransferF(sacZ, sacZamV, sacZphV, sacH2, sacH2amV, sacH2phV, Coh2, Adm2, Pha2);
	nfcorrV[1] = ncorr2;
	// and apply to the one with larger correlation coef first
	if( ncorr1 > ncorr2 ) {
		Coh2.clear(); Adm2.clear(); Pha2.clear();
		if(ncorr1>=nfcorr_min) { ApplyCorrection( sacZ, sacH1, Coh1, Adm1, Pha1 ); sacZamV.clear(); }
		ncorr1 = CalcTransferF(sacH2, sacH2amV, sacH2phV, sacH1, sacH1amV, sacH1phV, Coh1, Adm1, Pha1);
		nfcorrV[2] = ncorr1;
		if(ncorr1>=nfcorr_min) { ApplyCorrection( sacH2, sacH1, Coh1, Adm1, Pha1 ); sacH2amV.clear(); }
		ncorr1 = CalcTransferF(sacZ, sacZamV, sacZphV, sacH2, sacH2amV, sacH2phV, Coh1, Adm1, Pha1);
		nfcorrV[3] = ncorr1;
		if(ncorr1>=nfcorr_min) { ApplyCorrection( sacZ, sacH2, Coh1, Adm1, Pha1 ); sacZamV.clear(); }
	} else {
		Coh1.clear(); Adm1.clear(); Pha1.clear();
		if(ncorr2>=nfcorr_min) { ApplyCorrection( sacZ, sacH2, Coh2, Adm2, Pha2 ); sacZamV.clear(); }
		ncorr2 = CalcTransferF(sacH1, sacH1amV, sacH1phV, sacH2, sacH2amV, sacH2phV, Coh2, Adm2, Pha2);
		nfcorrV[2] = ncorr2;
		if(ncorr2>=nfcorr_min) { ApplyCorrection( sacH1, sacH2, Coh2, Adm2, Pha2 ); sacH1amV.clear(); }
		ncorr2 = CalcTransferF(sacZ, sacZamV, sacZphV, sacH1, sacH1amV, sacH1phV, Coh2, Adm2, Pha2);
		nfcorrV[3] = ncorr2;
		if(ncorr2>=nfcorr_min) { ApplyCorrection( sacZ, sacH1, Coh2, Adm2, Pha2 ); sacZamV.clear(); }
	}
	return nfcorrV;
}
*/

