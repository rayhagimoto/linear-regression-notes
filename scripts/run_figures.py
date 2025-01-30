from violations import HomoscedasticityStudy

def render(s):
    "Wrapper that runs simulation and renders a plot"
    study = s()
    study.simulate()
    study.render_plots()

if __name__ == '__main__':
    
    studies = [
        HomoscedasticityStudy,
    ]

    for study in studies:
        render(study)