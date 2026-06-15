{
	UInt_t messungen = 4;
	TString histogrammName = "Driftzeiten";
	TString dateien[4] = {"2500V.root", "2600V.root", "2700V.root", "2800V.root"};
	TString titel[4] = {"2500 V", "2600 V", "2700 V", "2800 V"};

	TLegend* leg = new TLegend(0.6, 0.5, 0.9, 0.7);
	leg->SetHeader("Spannungen");

	Bool_t first=true;
	UInt_t num = messungen;
	do {
		--num;
		TFile::Open(dateien[num]);
		TH1* plot = static_cast<TH1*>(gDirectory->FindObjectAny(histogrammName));
		plot->SetLineColor(num+1);
		if (first) {
			plot->Draw();
			first=false;
		} else {
			plot->Draw("same");
		}
		leg->AddEntry(plot, titel[num], "lep");
	} while (num != 0);

	leg->Draw("SAME");
}
