%#define OVTK_StimulationId_Target                            0x00008205 // 33285
%#define OVTK_StimulationId_NonTarget                         0x00008206 // 33286

function [box_out]=bci_Process(box_in)
global sync;
global events;


%out = myeeglogger(1,5,'RodrigoRamele',sprintf('p300_%f.dat',box_in.clock));
sync = sync + 1 ;

%box_in.inputs{1}.header

%box_in.inputs{1}.buffer{1}

disp( sprintf('Stimulation  Value: %d\n', sync) );

%for i = 1:size(box_in.inputs{1}.buffer,2)
    
%    chunk = box_in.inputs{1}.buffer{i};
    
    %chunk.start_time
    %chunk.end_time
    %chunk.matrix_data
    
%    chunk
%end


    % we iterate over the pending chunks on input 2 (STIMULATIONS)
    for i = 1: OV_getNbPendingInputChunk(box_in,2)

         % we pop the first chunk to be processed, note that box_in is used as
         % the output variable to continue processing
         [box_in, start_time, end_time, stim_set] = OV_popInputBuffer(box_in,2);

         
         % the stimulation set is a 3xN matrix.
         % The openvibe stimulation stream even sends empty stimulation sets
         % so the following boxes know there is no stimulation to expect in
         % the latest time range. These empty chunk are also in the matlab buffer.
         if(numel(stim_set) >= 3) % at least one stim in the set.
             fprintf('i:%d Start Time %f End_Time %f - Received stimulation code %i at time %f column/row %d \n', i,start_time, end_time, stim_set(1), stim_set(2), stim_set(3));
             %if stim_set(1) == box_in.settings(1).value % I'm lazy... I only check the first stim in the set.
             %    box_in.user_data.trigger_state = ~box_in.user_data.trigger_state;
             %    disp('Trigger is switched.')
             events = [ events [stim_set(1) ; stim_set(2)] ];
         end
    end
     
    
    for i = 1: OV_getNbPendingInputChunk(box_in,3)

         % we pop the first chunk to be processed, note that box_in is used as
         % the output variable to continue processing
         [box_in, start_time, end_time, stim_set] = OV_popInputBuffer(box_in,3);

         
         % the stimulation set is a 3xN matrix.
         % The openvibe stimulation stream even sends empty stimulation sets
         % so the following boxes know there is no stimulation to expect in
         % the latest time range. These empty chunk are also in the matlab buffer.
         if(numel(stim_set) >= 3) % at least one stim in the set.
             fprintf('i:%d Start Time %f End_Time %f - Received stimulation code %i at time %f \n', i,start_time, end_time, stim_set(1), stim_set(2));
             %if stim_set(1) == box_in.settings(1).value % I'm lazy... I only check the first stim in the set.
             %    box_in.user_data.trigger_state = ~box_in.user_data.trigger_state;
             %    disp('Trigger is switched.')
             %end
             events = [ events [stim_set(1) ; stim_set(2)] ];
         end
     end

box_out = box_in;

end

