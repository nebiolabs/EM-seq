def notificationHtml(String status) {
    def isSuccess = status == 'SUCCESS'
    def bannerBg = isSuccess ? '#dff0d8' : '#f2dede'
    def bannerBorder = isSuccess ? '#d6e9c6' : '#ebccd1'
    def bannerColor = isSuccess ? '#3c763d' : '#a94442'
    def bannerMsg = isSuccess ? 'Execution completed successfully!' : 'Execution failed!'

    def rows = []
    rows << ['Pipeline', params.workflow ?: '-']
    rows << ['Run name', workflow.runName]
    rows << ['Launch time', workflow.start]
    if (workflow.complete) {
        rows << ['Ending time', "${workflow.complete} (duration: ${workflow.duration})"]
    }
    def stats = null
    try { stats = workflow.stats } catch (Exception e) { }
    if (stats) {
        rows << ['Tasks stats', "Succeeded: ${stats.succeededCount} &nbsp; Cached: ${stats.cachedCount} &nbsp; Ignored: ${stats.ignoredCount} &nbsp; Failed: ${stats.failedCount}"]
    }
    rows << ['Launch directory', workflow.launchDir]
    rows << ['Work directory', workflow.workDir]
    rows << ['Project directory', workflow.projectDir]
    rows << ['Script name', workflow.scriptName]
    rows << ['Script ID', workflow.scriptId]
    rows << ['Workflow session', workflow.sessionId]
    if (workflow.profile) rows << ['Workflow profile', workflow.profile]
    rows << ['Nextflow version', "${workflow.nextflow.version}, build ${workflow.nextflow.build}"]

    def rowsHtml = rows.collect { pair ->
        "    <tr><td style=\"padding:4px 10px;color:#666;width:180px;vertical-align:top;\">${pair[0]}</td><td style=\"padding:4px 10px;word-break:break-all;\">${pair[1] ?: '-'}</td></tr>"
    }.join('\n')

    def errorBlock = ''
    if (!isSuccess && workflow.errorMessage) {
        errorBlock = """
  <p><b>Error:</b></p>
  <pre style="background:#f9f2f2;padding:10px;border-radius:4px;color:#a94442;white-space:pre-wrap;">${workflow.errorMessage}</pre>
"""
    }

    return """\
<!DOCTYPE html>
<html>
<head><meta charset="utf-8"></head>
<body style="font-family:Helvetica,Arial,sans-serif;max-width:800px;color:#333;padding:20px;">
  <h1 style="border-bottom:1px solid #ddd;padding-bottom:10px;">Workflow ${isSuccess ? 'completion' : 'failure'} notification</h1>
  <h2 style="margin-top:0;">Run Name: ${workflow.runName}</h2>

  <div style="background:${bannerBg};border:1px solid ${bannerBorder};color:${bannerColor};padding:10px 15px;border-radius:4px;margin:15px 0;">
    ${bannerMsg}
  </div>
${errorBlock}
  <p>The command used to launch the workflow was as follows:</p>
  <pre style="background:#f4f4f4;padding:10px;border-radius:4px;font-size:13px;white-space:pre-wrap;word-break:break-all;">${workflow.commandLine}</pre>

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

    def notified = false

    workflow.onError {
        def recipients = [params.email, params.admin_email].findAll { it }.join(',')
        if (!recipients) return
        sendMail(
            to: recipients,
            subject: "[Pipeline FAILED] ${params.workflow ?: 'unknown'} - ${workflow.runName}",
            body: notificationHtml('FAILED'),
            type: 'text/html'
        )
        notified = true
    }

    workflow.onComplete {
        if (notified) return
        def recipients = [params.email, params.admin_email].findAll { it }.join(',')
        if (!recipients) return
        def status = workflow.success ? 'SUCCESS' : 'FAILED'
        sendMail(
            to: recipients,
            subject: "[Pipeline ${status}] ${params.workflow} - ${workflow.runName}",
            body: notificationHtml(status),
            type: 'text/html'
        )
    }
}
