function [dataClass] = ApplyDataDictionary(dataPacket)
            
    dataLabels = {'Time_s'...
              'RED_Fx_Sat_N' ...
              'RED_Fy_Sat_N' ...
              'RED_Tz_Sat_Nm' ...
              'RED_Px_m' ...
              'RED_Py_m' ...
              'RED_Rz_rad' ...
              'RED_Vx_mpers' ...
              'RED_Vy_mpers' ...
              'RED_RzD_radpers' ...
              'RED_Ax_mpers2' ...
              'RED_Ay_mpers2' ...
              'RED_RzDD_radpers2' ...
              'BLACK_Fx_Sat_N' ...
              'BLACK_Fy_Sat_N' ...
              'BLACK_Tz_Sat_Nm' ...
              'BLACK_Px_m' ...
              'BLACK_Py_m' ...
              'BLACK_Rz_rad' ...
              'BLACK_Vx_mpers' ...
              'BLACK_Vy_mpers' ...
              'BLACK_RzD_radpers' ...
              'BLACK_Ax_mpers2' ...
              'BLACK_Ay_mpers2' ...
              'BLACK_RzDD_radpers2' ...
              'BLUE_Fx_Sat_N' ...
              'BLUE_Fy_Sat_N' ...
              'BLUE_Tz_Sat_Nm' ...
              'BLUE_Px_m' ...
              'BLUE_Py_m' ...
              'BLUE_Rz_rad' ...
              'BLUE_Vx_mpers' ...
              'BLUE_Vy_mpers' ...
              'BLUE_RzD_radpers' ...
              'ARM_Shoulder_Rz_rad' ...
              'ARM_Elbow_Rz_rad' ...
              'ARM_Wrist_Rz_rad' ...
              'ARM_Shoulder_RzD_radpers' ...
              'ARM_Elbow_RzD_radpers' ...
              'ARM_Wrist_RzD_radpers' ...
              'ARM_Command_Mode' ...
              'ARM_Shoulder_Command' ...
              'ARM_Elbow_Command'...
              'ARM_Wrist_Command'...
              'RED_Thruster_ON_Count'...
              'BLACK_Thruster_ON_Count'...
              'BLUE_Thruster_ON_Count'};
                          
        for sizePacket = 1:size(dataPacket,2)
            try
                eval(['dataClass.' dataLabels{sizePacket} ' = dataPacket(:,sizePacket);'])
            catch
                eval(['dataClass.CustomUserData' num2str(sizePacket) ' = dataPacket(:,sizePacket);'])
            end
        end
        
        
                    
end