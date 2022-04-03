---
name: Outreachy Code Project
about: Template for Outreachy Applicants interested in the Code project
title: 'Outreachy Code Project: <Applicant>'
labels: ''
assignees: ''

---

Applicant: <applicant github handle>

Welcome to the Ersilia Open Source Initiative. This issue will serve to track all your contributions for the project “Improve the Ersilia Model Hub, a FOSS platform offering pre-trained AI/ML models for research”.

Please tick the tasks as you complete them. To make a final application it is not required to have completed all tasks. This project requires knowledge of the Python programming language. The tasks are not ordered from more to less important, they are simply related to different skills. Start where you feel most comfortable.

---
### Initial steps
- [ ] Record your application for the project in the Outreachy website referencing this issue. Please make sure to select the right project on the website.
- [ ] Join the Slack channel to follow public communications.
- [ ] Comment under this issue explaining why you are interested in this project.
---
### Installation of the Ersilia Model Hub
- [ ] Install the `ersilia` library.
- [ ] Add a screenshot under this issue showing that you are able to run one model (for example, the chemprop-antibiotic model)
- [ ] Fetch at least 3 models from the Ersilia Model Hub. You can find these models with the `ersilia catalog` command. Add a screenshot of the *local* catalog (`ersilia catalog –local`)
---
### CLI
- [ ] Check if there are open issues related to the command line interface. Continue with the next tasks if they are open.
- [ ] Select one issue related to improving the CLI and request to be assigned to it.
- [ ] Link the #PR as a comment under this issue.
- [ ] Make any changes required in the PR and tick this box once it has been approved.
- [ ] Suggest at least one missing feature in the CLI (one sentence is enough, for example: “Add command to estimate memory usage of a particular model”).
---
### Python library
- [ ] Add a screenshot showing that you are able to run predictions using `ersilia` as a Python library (find more information  [here](https://ersilia.gitbook.io/ersilia-book/quick-start/antibiotic-activity-prediction)). Ideally, use a Jupyter notebook.
- [ ] Create a simple [Streamlit](https://streamlit.io/) app using the `ersilia` Python library. The app can have an input and an output box, and perhaps a few models to select. Add a screenshot of the app as seen in your browser.
- [ ] Write a docstring for the ErsiliaModel class. Use the Google Python Style guide. Paste the docstring as a comment below (do not use a PR).
---
### Scientific content
- [ ]  Check the models available in the Hub
- [ ] Select one model from the list and write a technical card (what is the model for, what input, which data was used to create it, what kind of ML algorithm uses…) for it
- [ ] Add your card as a comment to this issue
- [ ] Search the scientific literature and suggest 3 new models (comment in this issue) that would be relevant to incorporate in the Hub.
---
### Other
If you have interest in working on related topics, or have new suggestions, please do the following
- [ ] Add a comment in this issue with your new idea, tagging the mentor 
- [ ] Get feedback from the mentor and act accordingly
- [ ] Link in the comments any other PR you have contributed to.
---
### Community
- [ ] Look up two other projects and comment on their issues with feedback on one of their tasks
- [ ] If you have feedback from your peers, answer it in this issue.
---
### Final application
- [ ] I have answered all comments from mentors and contributors
- [ ] All PR or issues assigned to me are complete
- [ ] I have submitted my final application to the project
