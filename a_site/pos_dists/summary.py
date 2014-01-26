def collect_data(hr):
  f = open(str(hr)+'h_pos_dists/'+str(hr)+'h_pos_dist_summary.txt', 'w')
  f.write('Size\tOffset.0\t Offset.1\tOffset.2\n')
  for size in range(26, 37):
    g = open(str(hr)+'h_pos_dists/'+str(size)+'_fprint_pos.txt', 'r')
    g.readline()
    dist = []
    for i in range(3):
      line = g.readline().strip().split('\t')
      dist.append(line[1])
    f.write(str(size)+'\t'+'\t'.join(dist)+'\n')
    g.close()
  f.close()

if __name__ == '__main__':
  for hr in [0, 3, 6, 9, 12, 15]:
    collect_data(hr)
