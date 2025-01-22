import json
import subprocess
from datetime import datetime

import requests

# REPO_OWNER = 'ersilia-os'
# REPO_NAME = 'ersilia'
REPO_OWNER = "azycn"
REPO_NAME = "alice.github.io"
AUTH_TOKEN = None


def send_exception_issue(E: Exception, lastinput: str):
    subprocess.run(["gh", "auth", "login"])  # user login
    auth_out = subprocess.run(
        ["gh", "auth", "status", "-t"], capture_output=True
    )  # retrieve user's auth token

    # if stdout not there, parse stderr for Token
    if len(auth_out.stdout) == 0:
        # parse stderr for Token
        output = str(auth_out.stderr).split()
        found_token = False
        for i in output:
            if found_token:
                AUTH_TOKEN = i[0:-2]
                break
            if i == "Token:":
                found_token = True
        # if token not found, print error
        if not found_token:
            print("No token found.")
            print(str(auth_out.stderr)[2 : len(str(auth_out.stdout)) - 1])
            return
    else:
        # parse stdout
        output = str(auth_out.stdout).split()

        found_token = False
        for i in output:
            if found_token:
                AUTH_TOKEN = i[0:-2]
                break
            if i == "Token:":
                found_token = True
        print(AUTH_TOKEN)
        # if token not found, print error
        if not found_token and AUTH_TOKEN == None:
            print("No token found.")
            print(str(auth_out.stderr)[2 : len(str(auth_out.stdout)) - 1])
            return

    url = "https://api.github.com/repos/%s/%s/issues" % (REPO_OWNER, REPO_NAME)

    # Create an authenticated session to create the issue
    session = requests.Session()
    # Create our issue and headers to authenticate/format
    headers = {
        "Authorization": "token %s" % AUTH_TOKEN,
        "Accept": "application/vnd.github+json",
    }

    # Building the message for the body of the GitHub Issue
    dt_only = datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    datetime_str = "Date and time of issue: " + dt_only
    exception_str = "Exception message: " + str(E)
    lastinput_str = "Last Input: " + (
        lastinput if (lastinput != "") else "not provided."
    )
    # include log information here
    body_message = (
        datetime_str + "\n \n \n" + exception_str + "\n \n \n" + lastinput_str
    )

    issue = {
        "title": "New issue: " + dt_only,
        "body": body_message,
        "labels": ["Reported Issue"],
    }

    # Post the request to create an issue
    r = session.post(url, data=json.dumps(issue), headers=headers)

    if r.status_code == 201:
        print("Created issue!")
    else:
        print("Couldn't create issue.")
        print("Response:", r.status_code, r.content)


# previous email reporting:

# import smtplib, ssl
# from datetime import datetime

# smtp_server = "smtp.gmail.com"
# port = 587  # For starttls
# sender_email = "ersilia.errors@gmail.com"
# password = 'xjaluhzxtrnazpro' # insecure: the password will be visible from the code. Is this permissible?
#     # note - turn on 2 step verification, generate an "app password", and replace the password with this app pass

# sender_email = "ersilia.errors@gmail.com"
# receiver_email = "ersilia.errors@gmail.com"

# def send_exception_report_email(E:Exception, log_text:str = ""):
#     # get date/time info for message
#     datetime_str = datetime.now().strftime("%d/%m/%Y %H:%M:%S")

#     # Create a secure SSL context
#     context = ssl.create_default_context()

#     # Try to log in to server and send email
#     try:
#         server = smtplib.SMTP(smtp_server,port)
#         server.starttls(context=context) # Secure the connection
#         server.login(sender_email, password)

#         message = f"Subject: Exception {type(E)}: {datetime_str}\n\n"
#         message += f"{str(E)}\n\n"
#         message += log_text
#         server.sendmail(sender_email, receiver_email, message)

#     except Exception as e:
#         # Print any error messages to stdout
#         print(e)
#     finally:
#         server.quit()
