#!/usr/bin/octave

pkg load io;

args = argv ();

if rows (args) ~= 1
  disp ('usage: plot_csv.m <csv-file>')
  exit
end

filename = args{1};

data = csv2cell (filename);

headers = data (1, :);
values_cells = data (2:end, :);
values = cell2mat (values_cells);

# start_n = 15;
# end_n = 150;

methods_count = columns (headers);
ns = values (:, 1);

figure (1, 'visible', 'off');
hold on;
for i = 2:methods_count
  header = headers{i};
  err = values (:, i);

  figure (1, 'visible', 'off');
  plot (ns, err);
end

legend (headers{2:end});
xlabel ('N');
ylabel ('error');
h = figure (1);
waitfor (h);

