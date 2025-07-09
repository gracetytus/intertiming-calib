#include <TFile.h>
#include <TH1D.h>
#include <iostream>

int main() {
    TFile* f = TFile::Open("p1_out.root", "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Failed to open file\n";
        return 1;
    }

    TH1D* h = (TH1D*)f->Get("tdiff_1_2");
    if (!h) {
        std::cerr << "Histogram not found\n";
        return 2;
    }

    std::cout << "Entries: " << h->GetEntries() << std::endl;

    f->Close();
    delete f;
    return 0;
}

