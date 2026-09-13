// `wf` is the live WorkflowMetadata object captured by registerEmailNotifications
def notificationHtml(String status, String pipelineName, wf) {
    def isSuccess = status == 'SUCCESS'
    def bannerBg = isSuccess ? '#dff0d8' : '#f2dede'
    def bannerBorder = isSuccess ? '#d6e9c6' : '#ebccd1'
    def bannerColor = isSuccess ? '#3c763d' : '#a94442'
    def bannerMsg = isSuccess ? 'Execution completed successfully!' : 'Execution failed!'

    def rows = []
    rows << ['Pipeline', pipelineName ?: '-']
    rows << ['Run name', wf.runName]
    rows << ['Launch time', wf.start]
    if (wf.complete) {
        rows << ['Ending time', "${wf.complete} (duration: ${wf.duration})"]
    }
    def stats = null
    try { stats = wf.stats } catch (Exception e) { }
    if (stats) {
        rows << ['Tasks stats', "Succeeded: ${stats.succeededCount} &nbsp; Cached: ${stats.cachedCount} &nbsp; Ignored: ${stats.ignoredCount} &nbsp; Failed: ${stats.failedCount}"]
    }
    rows << ['Launch directory', wf.launchDir]
    rows << ['Work directory', wf.workDir]
    rows << ['Project directory', wf.projectDir]
    rows << ['Script name', wf.scriptName]
    rows << ['Script ID', wf.scriptId]
    rows << ['Workflow session', wf.sessionId]
    if (wf.profile) rows << ['Workflow profile', wf.profile]
    rows << ['Nextflow version', "${wf.nextflow.version}, build ${wf.nextflow.build}"]

    def rowsHtml = rows.collect { pair ->
        "    <tr><td style=\"padding:4px 10px;color:#666;width:180px;vertical-align:top;\">${pair[0]}</td><td style=\"padding:4px 10px;word-break:break-all;\">${pair[1] ?: '-'}</td></tr>"
    }.join('\n')

    def errorBlock = ''
    if (!isSuccess && wf.errorMessage) {
        errorBlock = """
  <p><b>Error:</b></p>
  <pre style="background:#f9f2f2;padding:10px;border-radius:4px;color:#a94442;white-space:pre-wrap;">${wf.errorMessage}</pre>
"""
    }

    return """\
<!DOCTYPE html>
<html>
<head><meta charset="utf-8"></head>
<body style="font-family:Helvetica,Arial,sans-serif;max-width:800px;color:#333;padding:20px;">
  <h1 style="border-bottom:1px solid #ddd;padding-bottom:10px;">Workflow ${isSuccess ? 'completion' : 'failure'} notification</h1>
  <h2 style="margin-top:0;">Run Name: ${wf.runName}</h2>

  <div style="background:${bannerBg};border:1px solid ${bannerBorder};color:${bannerColor};padding:10px 15px;border-radius:4px;margin:15px 0;">
    ${bannerMsg}
  </div>
${errorBlock}
  <p>The command used to launch the workflow was as follows:</p>
  <pre style="background:#f4f4f4;padding:10px;border-radius:4px;font-size:13px;white-space:pre-wrap;word-break:break-all;">${wf.commandLine}</pre>

  <h2 style="border-bottom:1px solid #ddd;padding-bottom:10px;margin-top:30px;">Execution summary</h2>
  <table style="border-collapse:collapse;font-size:14px;">
${rowsHtml}
  </table>
</body>
</html>
"""
}

def registerEmailNotifications() {
    if (params.dry_run || workflow.stubRun) return

    // workflow.onError/onComplete run detached from this binding, so bare `params.*`/`workflow.*`
    // throw NPEs there; capture into locals below so closure capture picks them up instead.
    def wf = workflow
    def recipients = [params.email, params.admin_email].findAll { it }.join(',')
    def pipelineName = params.workflow ?: 'unknown'

    def notified = false

    workflow.onError {
        if (!recipients) return
        sendMail(
            to: recipients,
            subject: "[Pipeline FAILED] ${pipelineName} - ${wf.runName}",
            body: notificationHtml('FAILED', pipelineName, wf),
            type: 'text/html'
        )
        notified = true
    }

    workflow.onComplete {
        if (notified) return
        if (!recipients) return
        def status = wf.success ? 'SUCCESS' : 'FAILED'
        sendMail(
            to: recipients,
            subject: "[Pipeline ${status}] ${pipelineName} - ${wf.runName}",
            body: notificationHtml(status, pipelineName, wf),
            type: 'text/html'
        )
    }
}
