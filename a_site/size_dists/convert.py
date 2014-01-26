for hr in [9]:
  f = open(str(hr) + 'h_foot_size_dist.txt', 'r')
  data = []
  for line in f:
    line = line.strip().split()
    data.append([line[1], line[3]])
  f.close()
  g = open(str(hr) + 'h_foot_size_dist.txt', 'w')
  g.write('#Size, Read Count\n')
  for line in data:
    g.write('\t'.join(line)+'\n')
  g.close()
