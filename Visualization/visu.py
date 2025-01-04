




def init_plot():
    fig = plt.figure(figsize=(12,8))
    plt.set_cmap('jet')
    ax = fig.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
    ax[0,0].set_xlim(0, lx_domain)
    ax[0,0].set_ylim(0, ly_domain)
    ax[0,0].set_xlabel('x')
    ax[0,0].set_ylabel('y')
    time_template = 'time = %.2fs'
    time_text = ax[0,0].text(0.5*lx_domain,1.1*ly_domain,'') 
    time_text.set_text(time_template%(0.))
    scat1 = ax[0,0].scatter(pos[:,0],pos[:,1],c=pres)
    div = make_axes_locatable(ax[0,0])
    cax1 = div.append_axes('right', '5%', '5%')

    scat2 = ax[0,1].scatter(pos[:,0],pos[:,1],c=rho)
    div = make_axes_locatable(ax[0,1])
    cax2 = div.append_axes('right', '5%', '5%')

    quiv = ax[1,0].quiver(pos[:,0],pos[:,1], vel[:,0], vel[:,1], color=['r','b','g'], scale=21)

    scat3 = ax[1,1].scatter(pos[:,0],pos[:,1],c=np.sqrt(vel[:,0]**2+vel[:,1]**2))
    div = make_axes_locatable(ax[1,1])
    cax3 = div.append_axes('right', '5%', '5%')
 
    return fig, cax1, cax2, cax3, scat1, scat2, scat3, quiv, time_text, time_template