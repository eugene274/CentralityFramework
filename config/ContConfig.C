Int_t GetRefMultiplicity (DataTreeEvent* event)
{
    Int_t RefMult = 0;
    
//     RefMult = event->GetNTracks();
    
    for (int i=0;i<event->GetNTracks();i++)
    {
        DataTreeTrack* track = event -> GetTrack(i);
        Bool_t cut = true;
        cut = cut && track->GetChiSq(0)/track->GetNDF(0) < 3  ;
        if (!cut)  continue;
        RefMult++;
    }
    
    return RefMult;
}


Bool_t isSelectedEvent(DataTreeEvent* event)
{
    return true;
}