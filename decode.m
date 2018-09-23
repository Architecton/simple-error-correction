% function [izhod, crc] = naloga3(vhod, n, k)

% Decodes vhod encoded with Hamming code H(n, k) and returns the result in izhod.
% Computes the CRC code of vhod using the CRC-8-CCITT standard and returns the result in crc.

function [izhod, crc] = decode(vhod, n, k)
  % Compute y and x dimensions of the H matrix
  H_y = log2(n + 1);
  
  % Compute matrix H
  H = [];
  for num = 1:n
    H = [H, dec2bin(num, H_y)' - '0'];
  endfor
  H = flipud(H);
  
  % The columns forming the identity matrix are indexed in powers of two [2^0..2^(H_y - 1)].
  I = eye(H_y);
  
  % Remove columns that form the identity matrix and add identity matrix to the end.
  del_columns = [];
  for val = 0:H_y - 1
    del_columns = [del_columns, pow2(val)];
  endfor
  H(:, del_columns) = [];
  H = [H, I];

  % Compute number of packets presented at the input.
  num_packets = length(vhod)/n;
  
  % Go over packets and apply error detection and correction.
  % Add data bits to result (izhod).
  izhod = [];
  for packet = 1:n:num_packets*n
    next_packet = vhod(packet:packet+n - 1);
    syndrome = mod(next_packet * H', 2);
    e = ismember(H', syndrome, 'rows');
    next_packet = xor(next_packet, e');
    izhod = [izhod, next_packet(1:k)];
  endfor
  
  % // CRC ////////////////////////////////////////
  % Compute CRC of input (vhod) using the CRC-8-CCITT standard
  
  % Initialize shift Register.
  shift_reg = zeros(1, 8);
  % Go over input.
  for ix_vhod = 1:length(vhod)
    % Perform XOR with next input value and last bit in the shift register.
    xor_res = xor(vhod(ix_vhod), shift_reg(end));
    % Initialize next shift register.
    shift_reg_next = shift_reg;
    % Add XOR result to first bit.
    shift_reg_next(1) = xor_res;
    % Perform XOR operations on following two bits with first and second bits in previous shift register state.
    shift_reg_next([2, 3]) = xor(shift_reg([1, 2]), [xor_res, xor_res]);
    % Shift bits in register.
    shift_reg_next(4:end) = shift_reg(3:end - 1);
    shift_reg = shift_reg_next;
  endfor
  
  % Convert value to hexadecimal.
  crc = fliplr(shift_reg);
  crc = num2str(crc);
  crc(isspace(crc)) = '';
  crc = bin2dec(crc);
  crc = dec2hex(crc);
endfunction
