function ETP = inso_forces(t_geol)
load 'insolation_and_temperature/latitude_inso_99_89.mat'
    if t_geol >=  -524.3 & t_geol <= -514.3
        ETP = zeros(40,1);
        for latitude_id = 1:1:40
            inso_at_lati = latitude_inso_99_89(:,latitude_id)';

            time_for_lati = -524.3:0.01:-514.3;
            key_past_time = min (time_for_lati( (time_for_lati - t_geol) >= 0 ) ) ;
            key_future_time = max (time_for_lati( (time_for_lati - t_geol) <= 0 ) ) ;
            %%%% find key insolation indexes and fractional contribution
            key_past_index = find( time_for_lati == key_past_time ) ;
            key_future_index = find( time_for_lati == key_future_time ) ;
            dist_to_past = abs( key_past_time - t_geol ) ;
            dist_to_future = abs( key_future_time - t_geol );
            %%%% fractional contribution of each keyframe
            if dist_to_past + dist_to_future == 0
                contribution_past = 1 ;
                contribution_future = 0 ;
            else
                contribution_past = dist_to_future / ( dist_to_past + dist_to_future ) ;
                contribution_future = dist_to_past / ( dist_to_past + dist_to_future ) ;
            end
    
            ETP(latitude_id,1) = contribution_past .* inso_at_lati(1,key_past_index) + contribution_future .* inso_at_lati(1,key_future_index);
            ETP(latitude_id,1) = ETP(latitude_id,1) - mean(inso_at_lati(1,:),2);
        end
    end