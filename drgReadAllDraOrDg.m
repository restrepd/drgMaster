function handles=drgReadAllDraOrDg(handles)
%VERY IMPORTANT: Read all the data from the dg file!!
%reading with an offset will not work

if handles.drg.dgordra==1
    %This is dra
    if exist(handles.drg.drta_p.fullName,'file')==2
        handles.data_dra=[];
        handles.data_dra=drgReadAllDra(handles.drg.drta_p.fullName,1,handles.drg);
    end
else
    %This is dg or rhd
    if exist(handles.drg.drta_p.fullName,'file')==2
        handles.data_dg=[];
        handles.data_dg=drgReadAllDg(handles.drg.drta_p.fullName,1,handles.drg);
    end
end