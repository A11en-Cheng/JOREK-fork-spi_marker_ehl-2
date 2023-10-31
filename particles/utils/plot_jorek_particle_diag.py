# ----------------------------------------------------------- #
# The function plots the particle diagnostics produce by a    #
# JOREK particle diagnostics                                  #
# ----------------------------------------------------------- #

# read a jorek particle diagnistics file. Zero ended bytes objects
# (S-type) are identified and transformed in strings
def read_jorek_particle_diagnostics_file(filename,filepath,separator):
  from h5py import File as h5pyfile
  from numpy import array,transpose
  groups = []
  fhandler = h5pyfile("".join([filepath,separator,filename]))
  for group in fhandler['groups'].values():
    diagnostics = {}
    for diag_name,diag in group.items():
      diagnostics[diag_name] = array(diag)
      if(len(diagnostics[diag_name].shape)>1):
        diagnostics[diag_name] = transpose(diagnostics[diag_name])
      if('S' in str(diagnostics[diag_name].dtype)):
        diagnostics[diag_name] = str(diagnostics[diag_name])
        diagnostics[diag_name] = diagnostics[diag_name][2:-1]
    groups.append(diagnostics)
  psi_axis = array(fhandler['psi_axis'])
  psi_bnd  = array(fhandler['psi_sep'])
  fhandler.close()
  return groups,psi_axis,psi_bnd

# plot profiles for multiple groups
def plot_time_profiles_multiple_groups(groups,xkeys,ykeys,titles,xlabels=[],ylabels=[],ncols=3,\
fontsize=18,linewidth=3,ploterror=False,addlegend=False):
  from numpy import ceil as npceil
  from numpy import abs as npabs
  from matplotlib.pyplot import subplots,show
  nykeys = len(ykeys)
  nplots = len(xkeys)*nykeys
  nrows = max(int(npceil(nplots/ncols)),1)
  if(nplots < ncols):
    ncols = nplots
  fig,axs = subplots(nrows=nrows,ncols=ncols,facecolor='white',edgecolor='white')
  for xkeyid,xkey in enumerate(xkeys):
    for ykeyid,ykey in enumerate(ykeys):
      for groupid,group in enumerate(groups):
        for particle_id,particle_prof in enumerate(group[ykey]):
          if(ploterror):
            if(particle_prof[0]!=0e0):
              axs[xkeyid*nykeys+ykeyid].plot(group[xkey],1e2*npabs((particle_prof-particle_prof[0])/particle_prof[0]),\
              linewidth=linewidth,label="".join(["group id: ",str(groupid),' particle id: ',str(particle_id)]))
          else:
            axs[xkeyid*nykeys+ykeyid].plot(group[xkey],particle_prof,linewidth=linewidth,\
            label="".join(["group id: ",str(groupid),' particle id: ',str(particle_id)]))
      if(ploterror):
        axs[xkeyid*nykeys+ykeyid].set_title("".join(["Error ",titles[xkeyid][ykeyid]]),fontsize=fontsize)
        axs[xkeyid*nykeys+ykeyid].set_ylabel("".join(["Error ",ylabels[ykeyid],"\%"]),fontsize=fontsize)
      else:
        axs[xkeyid*nykeys+ykeyid].set_title(titles[xkeyid][ykeyid],fontsize=fontsize)
        axs[xkeyid*nykeys+ykeyid].set_ylabel(ylabels[ykeyid],fontsize=fontsize)
      axs[xkeyid*nykeys+ykeyid].set_xlabel(xlabels[xkeyid],fontsize=fontsize)
      axs[xkeyid*nykeys+ykeyid].tick_params(axis='x',labelsize=fontsize)
      axs[xkeyid*nykeys+ykeyid].tick_params(axis='y',labelsize=fontsize)
      if(addlegend):
        axs[xkeyid*nykeys+ykeyid].legend(fontsize=fontsize)
      axs[xkeyid*nykeys+ykeyid].grid()
  show()

# main function
def read_analyse_plot_jorek_diagnostics(filename,filepath,separator,xkeys_prof,ykeys_prof,titles_prof,\
xlabels_prof,ylabels_prof,ncols,fontsize,linewidth,ploterror,addlegend):
  # read jorek particle diagnostics file
  p_groups,psi_axis,psi_bnd = read_jorek_particle_diagnostics_file(filename,filepath,separator)
  # ploty profiles for all groups
  plot_time_profiles_multiple_groups(p_groups,xkeys_prof,ykeys_prof,titles_prof,xlabels=xlabels_prof,\
  ylabels=ylabels_prof,ncols=ncols,fontsize=fontsize,linewidth=linewidth,ploterror=ploterror,\
  addlegend=addlegend)

# argument parser
def generate_argument_parser():
  from argparse import ArgumentParser
  parser = ArgumentParser(description='read and plot JOREK particle diagnostics')
  parser.add_argument('--filename','-f',type=str,action='store',required=False,\
  dest='filename',default='part_diag.h5',\
  help='name of the particle diagnostic file, default: part_diag.h5')
  parser.add_argument('--filepath','-fpath',type=str,action='store',required=False,\
  dest='filepath',default='.',help='path to the particle diagnostic file, default: .')
  parser.add_argument('--separator','-sep',type=str,action='store',required=False,\
  dest='separator',default='/',help='file separator, default: /')
  parser.add_argument('--xkeysprof','-xkp',type=str,nargs="*",action='store',required=False,\
  dest='xkeys_prof',default=['t'],help='coordinate values to plot, default: [t]')
  parser.add_argument('--ykeysprof','-ykp',type=str,nargs="*",action='store',required=False,\
  dest='ykeys_prof',default=['e','p_phi'],help='profiles to plot, default: [e,p_phi]')
  parser.add_argument('--titlesprof','-titp',type=str,nargs="*",action='store',required=False,\
  dest='titles_prof',default=[['Particle total energy vs time','Particle torroidal canonical momentum vs time']],\
  help='plot titles, default: [[Particle total energy vs time, Particle torroidal canonical momentum vs time]]')
  parser.add_argument('--xlabelsprof','-xlabp',type=str,nargs="*",action='store',required=False,\
  dest='xlabels_prof',default=['time [s]'],help='x labels titles, default: [time [s]]')
  parser.add_argument('--ylabelsprof','-ylabp',type=str,nargs="*",action='store',required=False,\
  dest='ylabels_prof',default=['Etot','P\phi'],help='y labels titles, default: [Etot [eV],P\phi]')
  parser.add_argument('--ncols',type=int,action='store',required=False,\
  dest='ncols',default=3,help='number of columns for subplot, default: 3')
  parser.add_argument('--fontsize',type=int,action='store',required=False,\
  dest='fontsize',default=18,help='plot font size, default: 18')
  parser.add_argument('--linewidth','-lwidth',type=int,action='store',required=False,\
  dest='linewidth',default=3,help='line width for subplot, default: 3')
  parser.add_argument('--ploerror','-plterr',type=bool,action='store',required=False,\
  dest='ploterror',default=True,help='plot profile error w.r.t. initial error, default: True')
  parser.add_argument('--addlegend','-al',type=bool,action='store',required=False,\
  dest='addlegend',default=False,help='add legend to plots, default: False')
  return parser.parse_args()

# Run main -------------------------------------------------- #
if __name__ == "__main__":
  args = generate_argument_parser()
  read_analyse_plot_jorek_diagnostics(args.filename,args.filepath,args.separator,\
  args.xkeys_prof,args.ykeys_prof,args.titles_prof,args.xlabels_prof,args.ylabels_prof,\
  args.ncols,args.fontsize,args.linewidth,args.ploterror,args.addlegend)
# ----------------------------------------------------------- #


