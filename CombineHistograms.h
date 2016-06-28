// Takes 2 1-D Histograms and a 2-D Histogram and creates a nice layout with
// The 2-D histogram surrounded by the 1-D Histograms
// This works by taking an image of the Y-Axis Histogram and then rotating it by 90˚
// Then displays it next to the 2-D histogram
//
// Call by running:
// CombineHistograms(<args>);
//
// Nota Bene: Arguments are separated by commas i.e.: CombineHistograms(arg1, arg2...)
// N.B.: Include this header file i.e.: #include "CombineHistograms.h"
//
// Required Arguments:
// X-Axis 1-D Histogram (TH1*)
// Y-Axis 1-D Histogramn (TH1*)
// 2-D Histogram (combination of those 2) (TH2*)
// Canvas To display it in (TCanvas*)
//
// Optional Arguments:
// Final Canvas Width (default: 1242) (int)
// Final Canvas Height (default: 755) (int)
// X-Axis Pad Description (default: X-Axis Pad) (const char)
// Y-Axis Pad Description (default: Y-Axis Pad) (const char)
//
// Grant Gollier

#include "TH2.h"
#include "TPad.h"
#include "TImage.h"

void CombineHistograms(TH1* HistoX, TH1* HistoY, TH2* Histo2D, TCanvas* CombinedCanvas,
                       //const char *XName = "PT", const char *YName = "Eta", const char *MainTitle = "Eta vs PT",
                       int Width = 1200, int Height = 740,
                       const char *XDescription = "X-Axis Pad", const char *YDescription = "Y-Axis Pad"){
    
    // Declare Canvases
    TCanvas *YCanvas = new TCanvas("YCanvas", "Y", 1000, 1000);
    //TCanvas *CombinedCanvas = new TCanvas("CombinedCanvas", "C");
    
    
    // Set Formating for the three histograms
    HistoY->GetXaxis()->SetTitle("");
    //HistoY->GetXaxis()->SetTitleSize(0.05);
    HistoY->GetYaxis()->CenterTitle();
    //HistoY->GetYaxis()->SetTitleSize(0.05);
    //HistoY->SetLabelSize(0.05, "xy");
    
    HistoX->GetXaxis()->SetTitle("");
    //HistoX->GetXaxis()->SetTitleSize(0.05);
    HistoX->GetYaxis()->CenterTitle();
    //HistoX->GetYaxis()->SetTitleSize(0.05);
    //HistoX->SetLabelSize(0.05, "xy");
    
    //Histo2D->SetXTitle(XName);
    //Histo2D->SetYTitle(YName);
    //Histo2D->SetTitle(MainTitle);
    Histo2D->GetXaxis()->CenterTitle();
    //Histo2D->GetXaxis()->SetTitleSize(0.05);
    Histo2D->GetYaxis()->CenterTitle();
    //Histo2D->GetYaxis()->SetTitleSize(0.05);
    //Histo2D->SetLabelSize(0.05, "xy");
    //Histo2D->SetStats(0);
    
    // Draw the Y-Axis Histogram so an image can be made from it
    YCanvas->cd();
    TPad *tempPad = new TPad("temppad", "temppad", 0.0, 0.0, 1.0, 1.0);
    tempPad->Draw();
    tempPad->cd();
    HistoY->UseCurrentStyle();
    HistoY->GetYaxis()->SetTitleOffset(2.2);
    HistoY->Draw();
    YCanvas->Update();
    tempPad->Update();
    
    // Draw the Canvas that will hold the 3 Histograms
    CombinedCanvas->Draw();
    CombinedCanvas->cd();
    CombinedCanvas->SetCanvasSize(Width, Height);
    
    // Pads for containng the 3 histograms
    TPad *XPad = new TPad("XPad", XDescription,0.4, 0.0, 1.0, 0.4);
    TPad *YPad = new TPad("YPad", YDescription,0.1028502,0.3258278,0.4522866,0.9748344);
    TPad *_2DPad = new TPad("2d", "2D plot pad",0.4, 0.4, 1.0, .9);
    
    // Draw the pads
    YPad->Draw();
    XPad->Draw();
    _2DPad->Draw();
    
    // Declare image
    TImage *img = TImage::Create();
    
    // Change to the respective pad and draw
    XPad->cd();
    HistoX->UseCurrentStyle();
    HistoX->Draw();
    
    YPad->cd();
    img->FromPad(tempPad);
    YCanvas->Close();
    img->Flip(90); // Rotates image by 90˚ CCW
    img->Draw();
    YPad->Update();
    
    _2DPad->cd();
    Histo2D->UseCurrentStyle();
    Histo2D->Draw("COLZ"); // Uses Color Plot with Key
}
