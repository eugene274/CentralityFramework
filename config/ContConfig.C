Int_t ContConfig(DataTreeEvent* event)
{
    Int_t RefMult;
    
    RefMult = event->GetNTracks();
    
    return RefMult;
}


Bool_t isSelectedEvent(DataTreeEvent* event)
{
    return true;
}